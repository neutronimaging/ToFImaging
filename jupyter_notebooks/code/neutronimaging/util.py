#!/usr/env/bash
# -*- coding: utf-8 -*-

"""
Provide various useful helper functions that
handles system level tasks.
"""
import os
from typing import Generator


def in_jupyter() -> bool:
    """check if current kernel is running as notebook backend"""
    try:
        from IPython import get_ipython

        kernel_name = get_ipython().__class__.__name__
        state = True if "ZMQ" in kernel_name else False
    except NameError:
        state = False
    return state


def probe_folder(root: str = ".") -> dict:
    """return folder structure as a dictionary"""
    return {
        os.path.basename(root): [
            os.path.join(root, me)
            if os.path.isfile(os.path.join(root, me))
            else probe_folder(os.path.join(root, me))
            for me in os.listdir(root)
        ]
    }


def _flatten_str_list(inlist: list) -> Generator:
    """Flatten a n-dimension nested list"""
    for item in inlist:
        if isinstance(item, str):
            yield item
        else:
            yield from _flatten_str_list(item)


def dir_tree_to_list(dir_tree: dict, flatten=True, sort=True) -> list:
    """Convert a dir tree (dict) to nested list"""
    _imglist = []
    for k, v in dir_tree.items():
        _imglist += [
            me if not isinstance(me, dict) else dir_tree_to_list(me) for me in v
        ]
    _imglist = list(_flatten_str_list(_imglist)) if flatten else _imglist
    return sorted(_imglist) if sort else _imglist


def convert_epics_timestamp_to_rfc3339_timestamp(epics_timestamp: float) -> float:
    # TIFF files from CG1D have EPICS timestamps.  From the Controls
    # Wiki:
    #
    # > EPICS timestamp. The timestamp is made when the image is read
    # > out from the camera. Format is seconds.nanoseconds since Jan 1st
    # > 00:00 1990.

    # Convert seconds since "EPICS epoch" to seconds since the "UNIX
    # epoch" so that Python can understand it.  I got the offset by
    # calculating the number of seconds between the two epochs at
    # https://www.epochconverter.com/
    EPOCH_OFFSET = 631152000
    unix_epoch_timestamp = EPOCH_OFFSET + epics_timestamp

    return unix_epoch_timestamp


if __name__ == "__main__":
    pass
