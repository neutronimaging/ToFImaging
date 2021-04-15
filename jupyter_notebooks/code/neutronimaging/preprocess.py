#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module contains necessary preprocessing toolkits for neutron imaging, including
- CG1D: reactor imaging beamline
-  
"""

import warnings
import json
import itertools
import pandas as pd
from PIL import Image
from datetime import datetime
from neutronimaging.npmath import find_edges_1d
from neutronimaging.util import dir_tree_to_list
from neutronimaging.util import probe_folder
from neutronimaging.util import convert_epics_timestamp_to_rfc3339_timestamp
from typing import List, Tuple
from typing import Union


def extract_metadata_tiff(tiffname: str) -> Tuple[list, list]:
    # default offset from the TIFF file
    # - str entry
    DATAACQMODE = 65018  # 'DataAcqModeStr:White Beam',
    DATATYPE = 65019  # 'DataTypeStr:Raw'
    DETECTOR_MANUFACTURER = 65026  # 'ManufacturerStr:Andor'
    FILENAMESTR = 65010  # 'FileNameStr:sam1_Cineole_LightOn'
    INSTRUMENT = 65011  # 'InstrumentStr:CG1D'
    MODELSTRING = 65025  # 'ModelStr:DW936_BV'

    # - float entry
    APERTURE_HR = 65068  # 'MotSlitHR.RBV:40.000000',
    APERTURE_HL = 65070  # 'MotSlitHL.RBV:40.000000'
    APERTURE_VT = 65066  # 'MotSlitVT.RBV:40.000000',
    APERTURE_VB = 65064  # 'MotSlitVB.RBV:40.000000',
    EXPOSURE_TIME = 65027  # 'ExposureTime:30.000000',
    IMGCOUNTER = 65031  # 'ImageCounter:77'
    TIME_SEC = 65002  # time secs
    TIME_NSEC = 65003  # time nano secs
    TIME_FULL = 65000  # full time
    IPTS = 65012  # 'IPTS:20267',
    ITEMS = 65013  # 'ITEMS:67144'
    GROUPID = 65020  # 'GroupID:113424'

    _metadata = dict(Image.open(tiffname).tag_v2)

    # time entry requires some special handling
    try:
        time_stamp = _metadata[TIME_SEC] + _metadata[TIME_NSEC] * 1e-9
    except:
        time_stamp = _metadata[TIME_FULL]
    finally:
        time_stamp = convert_epics_timestamp_to_rfc3339_timestamp(time_stamp)

    time_stamp_user_format = datetime.fromtimestamp(time_stamp).strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    header = [
        "filename",
        "time_stamp",
        "time_stamp_user_format",
        "data_acq_mode",
        "data_type",
        "detector_manufacturer",
        "filename_base",
        "instrument",
        "model_string",
        "aperture_HR",
        "aperture_HL",
        "aperture_VT",
        "aperture_VB",
        "exposure_time",
        "image_counter",
        "IPTS",
        "items",
        "groupid",
    ]

    data = [
        tiffname,
        time_stamp,
        time_stamp_user_format,
        _metadata[DATAACQMODE].split(":")[-1],
        _metadata[DATATYPE].split(":")[-1],
        _metadata[DETECTOR_MANUFACTURER].split(":")[-1],
        _metadata[FILENAMESTR].split(":")[-1],
        _metadata[INSTRUMENT].split(":")[-1],
        _metadata[MODELSTRING].split(":")[-1],
        float(_metadata[APERTURE_HR].split(":")[-1]),
        float(_metadata[APERTURE_HL].split(":")[-1]),
        float(_metadata[APERTURE_VT].split(":")[-1]),
        float(_metadata[APERTURE_VB].split(":")[-1]),
        float(_metadata[EXPOSURE_TIME].split(":")[-1]),
        int(_metadata[IMGCOUNTER].split(":")[-1]),
        int(_metadata[IPTS].split(":")[-1]),
        int(_metadata[ITEMS].split(":")[-1]),
        int(_metadata[GROUPID].split(":")[-1]),
    ]

    return data, header


def generate_config_CG1D(
    image_dir: Union[str, List],
    openbeam_dir: str,
    darkfield_dir: str,
    output: str = None,
    tolerance_aperature: float = 1.0,  # in mm
) -> dict:
    """frontend to allow list of rootdirs"""
    cfg_dict = {}
    if isinstance(image_dir, str):
        cfg_dict, df = _generate_config_CG1D(
            image_dir, openbeam_dir, darkfield_dir, None, tolerance_aperature
        )
    elif isinstance(image_dir, list):
        df_list = []
        for this_dir in image_dir:
            cfg_dict[this_dir], df = _generate_config_CG1D(
                this_dir, openbeam_dir, darkfield_dir, None, tolerance_aperature
            )
            df_list.append(df)
        df = pd.concat(df_list)
    else:
        raise ValueError(f"input dir has to be a string a list of strings")

    # dump dict to desired format if output file name provided
    if output is not None:
        _write_config_to_disk(cfg_dict, output, df)

    return cfg_dict


def _generate_config_CG1D(
    image_dir: str,
    openbeam_dir: str,
    darkfield_dir: str,
    output: str = None,
    tolerance_aperature: float = 1.0,  # in mm
) -> Tuple[dict, pd.DataFrame]:
    # build the metadata DataFrame
    img_list = []
    for _dir in (image_dir, openbeam_dir, darkfield_dir):
        img_list += [
            me
            for me in dir_tree_to_list(probe_folder(_dir), flatten=True, sort=True)
            if ".tif" in me.lower()
        ]
    meta_data = (extract_metadata_tiff(me) for me in img_list)

    # NOTE:
    # If the image list gets way too long, we can consider using multiprocessing
    # to speed up the parsing process as it is mostly an IO/thread bound process.
    md_data = [md for md, _ in meta_data]
    _, header = extract_metadata_tiff(img_list[0])
    df = pd.DataFrame(data=md_data, columns=header)

    # need to add a few extra feature vector for clustering
    lbs = ["aperture_HR", "aperture_HL", "aperture_VT", "aperture_VB"]
    lbs_binned = [f"{lb}_binned" for lb in lbs]
    for lb, lb_binned in zip(lbs, lbs_binned):
        df.loc[:, lb_binned] = 0

    # group by
    # - exposure_time
    # - (detector_name, aperture_[HR|HL|VT|VB])
    # and populate dict
    exposure_times = df["exposure_time"].unique()
    cfg_dict = {}
    for exposure in exposure_times:
        data_dict = cfg_dict[exposure] = {}
        # first, we need to bin the interval to form the category
        for lb, lb_binned in zip(lbs, lbs_binned):
            vals = df.loc[df["exposure_time"] == exposure, lb].unique()
            bin_edges = list(find_edges_1d(vals, atol=tolerance_aperature))
            for _low, _up in bin_edges:
                df.loc[
                    (df["exposure_time"] == exposure) & (df[lb].between(_low, _up)),
                    lb_binned,
                ] = df.loc[
                    (df["exposure_time"] == exposure) & (df[lb].between(_low, _up)),
                    lb,
                ].mean()
        # second, find the categories
        detector_names = df.loc[
            df["exposure_time"] == exposure, "detector_manufacturer"
        ].unique()
        aperture_HRs = df.loc[
            df["exposure_time"] == exposure, "aperture_HR_binned"
        ].unique()
        aperture_HLs = df.loc[
            df["exposure_time"] == exposure, "aperture_HL_binned"
        ].unique()
        aperture_VTs = df.loc[
            df["exposure_time"] == exposure, "aperture_VT_binned"
        ].unique()
        aperture_VBs = df.loc[
            df["exposure_time"] == exposure, "aperture_VB_binned"
        ].unique()
        categories = itertools.product(
            detector_names, aperture_HRs, aperture_HLs, aperture_VTs, aperture_VBs
        )
        # last, populate each categroy
        metadata_info_keys = [
            "detector_manufacturer",
            "aperture_HR",
            "aperture_HL",
            "aperture_VT",
            "aperture_VB",
        ]
        list_sample_keys = ["filename", "time_stamp", "time_stamp_user_format"]
        for i, me in enumerate(categories):
            _tmp = data_dict[f"config{i}"] = {}
            # generate metadata_infos (the common core)
            _tmp["metadata_infos"] = {k: v for k, v in zip(metadata_info_keys, me)}
            # generate list of images (data_type: Raw)
            # generate list of ob (data_type: OB)
            # generate list of df (data_type: DF)
            for groupname, datatype in zip(
                ("list_sample", "list_ob", "list_df"), ("Raw", "OB", "DF")
            ):
                _df_tmp = df.loc[
                    (df["exposure_time"] == exposure)
                    & (df["detector_manufacturer"] == me[0])
                    & (df["aperture_HR_binned"] == me[1])
                    & (df["aperture_HL_binned"] == me[2])
                    & (df["aperture_VT_binned"] == me[3])
                    & (df["aperture_VB_binned"] == me[4])
                    & (df["data_type"] == datatype),
                    list_sample_keys,
                ]
                _tmp[groupname] = [
                    {k: v for k, v in zip(list_sample_keys, row)}
                    for _, row in enumerate(_df_tmp.to_numpy())
                ]

            # generate for first images
            # generate for last images
            _tmp["first_images"] = {}
            _tmp["last_images"] = {}
            for groupname, datatype in zip(("sample", "ob", "df"), ("Raw", "OB", "DF")):
                _df_tmp = df.loc[
                    (df["exposure_time"] == exposure)
                    & (df["detector_manufacturer"] == me[0])
                    & (df["aperture_HR_binned"] == me[1])
                    & (df["aperture_HL_binned"] == me[2])
                    & (df["aperture_VT_binned"] == me[3])
                    & (df["aperture_VB_binned"] == me[4])
                    & (df["data_type"] == datatype),
                    list_sample_keys,
                ]
                _tmp["first_images"][groupname] = (
                    {
                        k: v
                        for k, v in zip(
                            list_sample_keys,
                            _df_tmp.sort_values("time_stamp").to_numpy()[0],
                        )
                    }
                    if _df_tmp.size > 0
                    else {}
                )
                _tmp["last_images"][groupname] = (
                    {
                        k: v
                        for k, v in zip(
                            list_sample_keys,
                            _df_tmp.sort_values("time_stamp").to_numpy()[-1],
                        )
                    }
                    if _df_tmp.size > 0
                    else {}
                )

            # generate time range
            # NOTE: need confirmation from original developer
            _tmp["time_range_s"] = {}
            first_sample_time = (
                _tmp["first_images"]["sample"]["time_stamp"]
                if "time_stamp" in _tmp["first_images"]["sample"].keys()
                else 0.0
            )
            first_ob_time = (
                _tmp["first_images"]["ob"]["time_stamp"]
                if "time_stamp" in _tmp["first_images"]["ob"].keys()
                else 0.0
            )
            _tmp["time_range_s"]["before"] = max(first_sample_time - first_ob_time, 0)
            last_sample_time = (
                _tmp["last_images"]["sample"]["time_stamp"]
                if "time_stamp" in _tmp["last_images"]["sample"].keys()
                else 0.0
            )
            last_ob_time = (
                _tmp["first_images"]["ob"]["time_stamp"]
                if "time_stamp" in _tmp["first_images"]["ob"].keys()
                else 0.0
            )
            _tmp["time_range_s"]["after"] = max(last_sample_time - last_ob_time, 0)

    # dump dict to desired format if output file name provided
    if output is not None:
        _write_config_to_disk(cfg_dict, output, df)

    return cfg_dict, df


def _write_config_to_disk(
    cfg_dict: dict,
    filename: str,
    dataframe: pd.DataFrame,
) -> None:
    _file_extension = filename.split(".")[-1]
    if "json" in _file_extension.lower():
        with open(filename, "w") as outputf:
            json.dump(cfg_dict, outputf, indent=2, sort_keys=True)
    elif "csv" in _file_extension.lower():
        dataframe.to_csv(filename, sep="\t", index=False)
    else:
        warnings.warn(
            f"Unsupported file extension provided: {_file_extension}, falling back to json"
        )
        filename += ".json"
        with open(filename, "w") as outputf:
            json.dump(cfg_dict, outputf, indent=2, sort_keys=True)


if __name__ == "__main__":
    pass
