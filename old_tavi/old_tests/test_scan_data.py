# -*- coding: utf-8 -*-
import copy

import numpy as np
import pytest

from tavi.data.scan_data import ScanData1D


@pytest.fixture
def scans1d():
    scan0001 = ScanData1D(x=np.array([0, 1, 2]), y=np.array([1, 2, 3]), norm=np.array([2, 2, 3]))
    scan0002 = ScanData1D(
        x=np.array([15, 15.1, 15.2, 15.3, 15.4, 15.1, 15.2, 15.3, 15.4, 15.5]),
        y=np.array([10, 12, 15, 42, 90, 31, 34, 105, 230, 3]),
        norm=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
    )
    scan0003 = ScanData1D(x=np.array([0.1, 1.1, 2.1]), y=np.array([1, 1, 1]))
    scan0004 = ScanData1D(x=np.array([-0.9, 0.1, 1.1, 2.1, 3.1]), y=np.array([10, 1, 1, 1, 10]))

    scans = (scan0001, scan0002, scan0003, scan0004)

    return scans


def test_scan_data_1d_renorm(scans1d):
    scan0001, *_ = scans1d
    scan0001.renorm(norm_val=2)
    assert np.allclose(scan0001.y, [1, 2, 2], atol=1e-3)
    assert np.allclose(scan0001.err, [1, np.sqrt(2), np.sqrt(3) / 3 * 2], atol=1e-3)


def test_rebin_grid_boundaries(scans1d):
    _, scan0002, *_ = scans1d

    test1 = copy.deepcopy(scan0002)

    test1.rebin_grid(rebin_params=(14.9, 15.6, 0.1))
    assert np.allclose(test1.x, np.array([14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6]), atol=1e-3)
    y_exp = np.array([np.nan, 10, (12 + 31) / 2, (15 + 34) / 2, (42 + 105) / 2, (90 + 230) / 2, 3, np.nan])
    assert np.allclose(test1.y, y_exp, atol=1e-3, equal_nan=True)
    err_exp = np.array(
        [
            np.nan,
            np.sqrt(10),
            np.sqrt(12 + 31) / 2,
            np.sqrt(15 + 34) / 2,
            np.sqrt(42 + 105) / 2,
            np.sqrt(90 + 230) / 2,
            np.sqrt(3),
            np.nan,
        ]
    )
    assert np.allclose(test1.err, err_exp, atol=1e-3, equal_nan=True)

    test2 = copy.deepcopy(scan0002)
    test2.rebin_grid(rebin_params=(15.1, 15.4, 0.1))
    assert np.allclose(test2.x, np.array([15.1, 15.2, 15.3, 15.4]), atol=1e-3)
    y_exp = np.array([(12 + 31) / 2, (15 + 34) / 2, (42 + 105) / 2, (90 + 230) / 2])
    assert np.allclose(test2.y, y_exp, atol=1e-3, equal_nan=True)
    err_exp = np.array(
        [
            np.sqrt(12 + 31) / 2,
            np.sqrt(15 + 34) / 2,
            np.sqrt(42 + 105) / 2,
            np.sqrt(90 + 230) / 2,
        ]
    )
    assert np.allclose(test2.err, err_exp, atol=1e-3, equal_nan=True)


def test_rebin_grid(scans1d):
    _, scan0002, *_ = scans1d
    scan0002.rebin_grid_renorm(
        rebin_params=(15.0, 15.5, 0.2),
        norm_col=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
        norm_val=4.0,
    )
    assert np.allclose(
        scan0002.x,
        [15.0, 15.2, 15.4],
        atol=1e-3,
    )
    assert np.allclose(
        scan0002.y,
        [
            10 / 2 * 4,
            (12 + 15 + 31 + 34) / (2 + 2 + 5 + 5) * 4,
            (42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
        ],
        atol=1e-3,
        equal_nan=True,
    )
    assert np.allclose(
        scan0002.err,
        [
            np.sqrt(10) / 2 * 4,
            np.sqrt(12 + 15 + 31 + 34) / (2 + 2 + 5 + 5) * 4,
            np.sqrt(42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
        ],
        atol=1e-3,
        equal_nan=True,
    )


def test_rebin_tol(scans1d):
    _, scan0002, *_ = scans1d

    scan0002.rebin_tol(
        rebin_params=(15.0, 15.5, 0.2),
        # weight_col=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
    )
    assert np.allclose(
        scan0002.x,
        [
            15,
            (2 * 15.1 + 5 * 15.1 + 2 * 15.2 + 5 * 15.2) / (2 + 2 + 5 + 5),
            (2 * 15.3 + 5 * 15.3 + 2 * 15.4 + 5 * 15.4) / (2 + 2 + 5 + 5),
        ],
        atol=1e-3,
    )
    assert np.allclose(
        scan0002.y,
        [10, (12 + 15 + 31 + 34) / 4, (42 + 90 + 105 + 230) / 4],
        atol=1e-3,
        equal_nan=True,
    )
    assert np.allclose(
        scan0002.err,
        [
            np.sqrt(10),
            np.sqrt(12 + 15 + 31 + 34) / 4,
            np.sqrt(42 + 90 + 105 + 230) / 4,
        ],
        atol=1e-3,
        equal_nan=True,
    )


def test_rebin_tol_renorm(scans1d):
    _, scan0002, *_ = scans1d

    scan0002.rebin_tol_renorm(
        rebin_params=(15.0, 15.5, 0.2),
        # norm_col=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
        norm_val=4.0,
    )
    assert np.allclose(
        scan0002.x,
        [
            15,
            (2 * 15.1 + 5 * 15.1 + 2 * 15.2 + 5 * 15.2) / (2 + 2 + 5 + 5),
            (2 * 15.3 + 5 * 15.3 + 2 * 15.4 + 5 * 15.4) / (2 + 2 + 5 + 5),
        ],
        atol=1e-3,
    )

    assert np.allclose(
        scan0002.y,
        [
            10 / 2 * 4,
            (12 + 15 + 31 + 34) / (2 + 2 + 5 + 5) * 4,
            (42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
        ],
        atol=1e-3,
        equal_nan=True,
    )
    assert np.allclose(
        scan0002.err,
        [
            np.sqrt(10) / 2 * 4,
            np.sqrt(12 + 15 + 31 + 34) / (2 + 2 + 5 + 5) * 4,
            np.sqrt(42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
        ],
        atol=1e-3,
        equal_nan=True,
    )


def test_sub(scans1d):
    scan0001, _, scan0003, *_ = scans1d
    scan_sub = scan0001 - scan0003
    assert np.allclose(scan_sub.x, [0, 1, 2])
    assert np.allclose(scan_sub.y, [0, 1, 2])
    assert np.allclose(scan_sub.err, [np.sqrt(2), np.sqrt(3), 2])


# def test_add(scans1d):
#     scan0001, _, scan0003, *_ = scans1d
#     scan_add = scan0001 + scan0003
#     assert np.allclose(scan_add.x, [0, 1, 2])
#     assert np.allclose(scan_add.y, [2, 3, 4])
#     assert np.allclose(scan_add.err, [np.sqrt(2), np.sqrt(3), 2])


# def test_sub_mismatch_x(scans1d):
#     scan0001, _, _, scan0004, *_ = scans1d
#     scan_sub = scan0001 - scan0004
#     assert np.allclose(scan_sub.x, [0, 1, 2])
#     assert np.allclose(scan_sub.y, [0, 1, 2])
#     assert np.allclose(scan_sub.err, [np.sqrt(2), np.sqrt(3), 2])
