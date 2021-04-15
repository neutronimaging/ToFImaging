#!/usr/env/bash
# -*- coding: utf-8 -*-

"""
This module contains various useful numpy-based math
subroutines, including
- find_edges_1d
"""

from typing import List
import numpy as np


def find_edges_1d(array: np.ndarray, atol: float = 1) -> List:
    """
    Locate the lower and upper edges of binning using given absoute
    tolerance.

    >> test_array = np.array([40, 40.5, 42, 41.7, 38])
    >> list(find_edges(test_array))
    >> [(37.5, 38.5), (38.5, 41.0), (41.0, 42.5)]
    """
    if len(array) <= 1:
        return [(min(array) - atol / 2, max(array) + atol / 2)]
    else:
        array = np.sort(array)
        gaps = array[1:] - array[:-1]
        lower_edge = [array[0] - atol / 2] + list(
            array[np.where(gaps > atol / 2)] + atol / 2
        )
        upper_edge = list(array[np.where(gaps > atol / 2)] + atol / 2) + [
            array[-1] + atol / 2
        ]
        return zip(lower_edge, upper_edge)


if __name__ == "__main__":
    pass
