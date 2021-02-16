#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
"""

import numpy as np
from typing import Type
from numpy.core.fromnumeric import sort
import pandas as pd
from NeuNorm.normalization import Normalization


def read_shutter_count(filename: str) -> pd.DataFrame:
    """Parse in shutter count data from csv"""
    _df = pd.read_csv(filename, sep="\t", names=["shutter_index", "shutter_counts"])
    _df = _df[_df["shutter_counts"] > 0]
    _df["shutter_n_ratio"] = _df["shutter_counts"] / _df["shutter_counts"].values[0]
    return _df


def read_shutter_time(filename: str) -> pd.DataFrame:
    """Parse in shutter time from csv"""
    _df = pd.read_csv(
        filename,
        sep="\t",
        names=["shutter_index", "start_frame", "end_frame"],
    )
    _df = _df[_df["end_frame"] > 0]
    # NOTE: the start/end frame here is delta, we need the absolute
    #       therefore, cumulative sum for times
    _lbs = ["start_frame", "end_frame"]
    _tmp = _df[_lbs].values
    _df[_lbs] = _tmp.flatten().cumsum().reshape(_tmp.shape)
    return _df


def read_spectra(filename: str) -> pd.DataFrame:
    """Parse in spectra data from csv"""
    return pd.read_csv(filename, sep="\t", names=["shutter_time", "counts"])


def merge_meta_data(
    shutter_count: pd.DataFrame,
    shutter_time: pd.DataFrame,
    spectra: pd.DataFrame,
) -> pd.DataFrame:
    """Consolidate meta data from three different dataframes into one"""
    _df = spectra.copy(deep=True)
    _df_shutter = pd.concat([shutter_count, shutter_time], axis=1)
    _df["run_num"] = spectra.index
    # initialize fields
    _df["shutter_index"] = -1
    _df["shutter_counts"] = -1
    for _, row in _df_shutter.iterrows():
        _idx, _cnt, _snr, _, _start, _end = row
        _df.loc[_df["shutter_time"].between(_start, _end), "shutter_index"] = int(_idx)
        _df.loc[_df["shutter_time"].between(_start, _end), "shutter_counts"] = int(_cnt)
        _df.loc[_df["shutter_time"].between(_start, _end), "shutter_n_ratio"] = _snr
    return _df


def skipping_meta_data(meta: pd.DataFrame) -> pd.DataFrame:
    """Skips first and last or each run in metadata"""
    _by_shutter = meta.groupby(['shutter_index'])
    # groupby returns (group_num, dataframe), thus the [1] first
    _with_skips = [item[1][1:-1] for item in _by_shutter]
    return pd.concat(_with_skips)


def load_images(raw_imamge_dir: str) -> Type[Normalization]:
    """Loading all Images into memory"""
    import glob
    from neutronimaging.util import in_jupyter

    o_norm = Normalization()

    # gather all image
    _img_names = [
        me for me in glob.glob(f"{raw_imamge_dir}/*.fits") if "_SummedImg" not in me
    ]
    _img_names.sort()

    o_norm.load(file=_img_names, notebook=in_jupyter())
    return o_norm


def calc_pixel_occupancy_probability(
    o_norm: Type[Normalization],
    metadata: pd.DataFrame,
) -> np.ndarray:
    """calculate pixel occupancy probability"""
    _imgs = np.array(o_norm.data["sample"]["data"])
    _pops = np.zeros_like(_imgs)

    # calculation is done on a per shutter index base
    for _idx in metadata["shutter_index"].unique():
        _run_num = metadata.loc[metadata["shutter_index"] == _idx, "run_num"].values
        _cnts = metadata.loc[metadata["shutter_index"] == _idx, "shutter_counts"].values
        _tmp = _imgs[_run_num, :, :].cumsum(axis=0)
        _pops[_run_num, :, :] = np.divide(_tmp, _cnts[:, np.newaxis, np.newaxis])
    return _pops


def correct_images(
    o_norm: Type[Normalization],
    metadata: pd.DataFrame,
    skip_first_and_last=False,
) -> np.ndarray:
    """
    Correct raw images based on shutter info in metadata
    """
    _img = np.array(o_norm.data["sample"]["data"])
    _pop = calc_pixel_occupancy_probability(o_norm, metadata)
    _snr = metadata["shutter_n_ratio"].values[:, np.newaxis, np.newaxis]
    _rst = _img / (1 - _pop) / _snr
    # NOTE: The very first and last image of each frame (shutter_index)
    #       needs specicial correction, therefore removing them from
    #       standard pipeline if specified
    if skip_first_and_last:
        _tmp = []
        for _idx in metadata["shutter_index"].unique():
            _run_num = metadata.loc[metadata["shutter_index"] == _idx, "run_num"].values
            _tmp += list(_run_num[1:-1])
        _idx_to_keep = np.array(_tmp)
        _rst = _rst[_idx_to_keep, :, :]
    return _rst


if __name__ == "__main__":
    import os

    _file_root = os.path.dirname(os.path.abspath(__file__))
    test_data_dir = os.path.join(_file_root, "../data")
    #
    shutter_counts_file = os.path.join(test_data_dir, "OB_1_005_ShutterCount.txt")
    df_shutter_count = read_shutter_count(shutter_counts_file)
    print(df_shutter_count)
    #
    shutter_time_file = os.path.join(test_data_dir, "OB_1_005_ShutterTimes.txt")
    df_shutter_time = read_shutter_time(shutter_time_file)
    print(df_shutter_time)
    #
    spectra_file = os.path.join(test_data_dir, "OB_1_005_Spectra.txt")
    df_spectra = read_spectra(spectra_file)
    print(df_spectra)
    #
    df_meta = merge_meta_data(df_shutter_count, df_shutter_time, df_spectra)
    print(df_meta)

    # test load images
    img_dir = test_data_dir
    o_norm = load_images(img_dir)
    print(type(o_norm))

    # test calculate pixel occupancy probability
    pop = calc_pixel_occupancy_probability(o_norm, df_meta)

    # test image correction
    imgs = correct_images(o_norm, df_meta)
