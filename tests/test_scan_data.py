# -*- coding: utf-8 -*-
import numpy as np
import pytest

from tavi.data.scan_data import ScanData1D


@pytest.fixture
def scans():
    scan0001 = ScanData1D(x=np.array([0, 1, 2]), y=np.array([1, 2, 3]))
    scan0002 = ScanData1D(
        x=np.array([15, 15.1, 15.2, 15.3, 15.4, 15.1, 15.2, 15.3, 15.4, 15.5]),
        y=np.array([10, 12, 15, 42, 90, 31, 34, 105, 230, 3]),
    )

    return (
        scan0001,
        scan0002,
    )


def test_scan_data_1d_renorm(scans):
    scan0001, *_ = scans
    scan0001.renorm(norm_col=np.array([2, 2, 3]), norm_val=2)
    assert np.allclose(scan0001.y, [1, 2, 2], atol=1e-3)
    assert np.allclose(scan0001.err, [1, np.sqrt(2), np.sqrt(3) / 3 * 2], atol=1e-3)


def test_rebin_grid_renorm(scans):
    _, scan0002, *_ = scans

    scan0002.rebin_grid_renorm(
        rebin_params=(15.0, 15.5, 0.2),
        norm_col=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
        norm_val=4.0,
    )
    assert np.allclose(
        scan0002.x,
        [15.1, 15.3, 15.5],
        atol=1e-3,
    )
    assert np.allclose(
        scan0002.y,
        [
            (10 + 12 + 15 + 31 + 34) / (2 + 2 + 2 + 5 + 5) * 4,
            (42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
            (3) / (5) * 4,
        ],
        atol=1e-3,
    )
    assert np.allclose(
        scan0002.err,
        [
            np.sqrt(10 + 12 + 15 + 31 + 34) / (2 + 2 + 2 + 5 + 5) * 4,
            np.sqrt(42 + 90 + 105 + 230) / (2 + 5 + 2 + 5) * 4,
            np.sqrt(3) / (5) * 4,
        ],
        atol=1e-3,
    )


def test_rebin_tol_renorm(scans):
    _, scan0002, *_ = scans

    scan0002.rebin_tol_renorm(
        rebin_params=(15.0, 15.5, 0.2),
        norm_col=np.array([2, 2, 2, 2, 2, 5, 5, 5, 5, 5]),
        norm_val=4.0,
    )
    assert np.allclose(
        scan0002.x,
        [
            (2 * 15 + 2 * 15.1 + 5 * 15.1 + 2 * 15.2 + 5 * 15.2) / (2 + 2 + 5 + 2 + 5),
            (2 * 15.3 + 5 * 15.3 + 2 * 15.4 + 5 * 15.4) / (2 + 2 + 5 + 5),
            (5 * 15.5) / (5),
        ],
        atol=1e-3,
    )
