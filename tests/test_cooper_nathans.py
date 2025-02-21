import numpy as np
import pytest

from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.sample import Sample

np.set_printoptions(floatmode="fixed", precision=4)


def test_copper_nathans_localQ(tas_params):
    tas, ei, ef, hkl, _, r0 = tas_params

    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, projection=None, R0=r0)
    mat = np.array(
        [
            [9583.2881, -4671.0614, -0.0000, 986.5610],
            [-4671.0614, 21359.2992, 0.0000, -4129.1553],
            [0.0000, 0.0000, 77.7036, 0.0000],
            [986.5610, -4129.1553, -0.0000, 864.3494],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_copper_nathans_hkl(tas_params):
    tas, ei, ef, hkl, _, r0 = tas_params

    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, R0=r0)
    mat = np.array(
        [
            [33305.0843, 33224.4963, -2651.8290, -5152.9962],
            [33224.4963, 33305.2609, -2651.8526, -5153.0102],
            [-2651.8290, -2651.8526, 1983.2037, 448.8024],
            [-5152.9962, -5153.0102, 448.8024, 864.3494],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_copper_nathans_projection(tas_params):
    tas, ei, ef, hkl, projection, r0 = tas_params

    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, projection=projection, R0=r0)
    mat = np.array(
        [
            [1.3306e05, -5.3037e03, -1.7660e-01, -1.0306e04],
            [-5.3037e03, 1.9832e03, 2.3558e-02, 4.4880e02],
            [-1.7660e-01, 2.3558e-02, 1.6135e02, 1.4003e-02],
            [-1.0306e04, 4.4880e02, 1.4003e-02, 8.6435e02],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_copper_nathans_list(tas_params):
    tas, ei, ef, _, _, r0 = tas_params

    rez_list = tas.rez(hkl_list=[(0, 0, 3), (0, 0, -3)], ei=[ei, ei + 1], ef=ef, projection=None, R0=r0)
    mat = np.array(
        [
            [9583.2881, -4671.0614, -0.0000, 986.5610],
            [-4671.0614, 21359.2992, 0.0000, -4129.1553],
            [0.0000, 0.0000, 77.7036, 0.0000],
            [986.5610, -4129.1553, -0.0000, 864.3494],
        ]
    )
    assert len(rez_list) == 4
    assert np.allclose(rez_list[0].mat, mat, atol=1e-1)
    assert np.allclose(rez_list[1].mat, mat, atol=1e-1)


@pytest.fixture
def tas_params():
    # cooper_nathans_CTAX

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = CooperNathans(fixed_ef=4.8)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Sample.from_json(sample_json_path)
    tas.mount_sample(sample)

    ei = 4.8
    ef = 4.8
    hkl = (0, 0, 3)

    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = False

    tas_params = (tas, ei, ef, hkl, projection, R0)

    return tas_params
