import numpy as np
import pytest

from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.instrument.tas import TAS
from tavi.sample import Sample


def test_validate_instrument_parameter(ctax):
    cn = CooperNathans(instrument=ctax)
    cn.validate_instrument_parameters()


def test_generate_hkle_list(ctax):
    cn = CooperNathans(instrument=ctax)
    hkle_list = cn.generate_hkle_list(
        hkl=[(h, h, 3) for h in np.arange(0, 1, 0.1)],
        en=[en for en in np.arange(0, 5, 0.5)],
    )
    assert len(hkle_list) == 100


def test_local_q(ctax):
    hkl = (0, 0, 3)
    R0 = False
    rez = ctax.cooper_nathans(hkl=hkl, en=0, projection=None, R0=R0)
    mat = np.array(
        [
            [9583.2881, -4671.0614, -0.0000, 986.5610],
            [-4671.0614, 21359.2992, 0.0000, -4129.1553],
            [0.0000, 0.0000, 77.7036, 0.0000],
            [986.5610, -4129.1553, -0.0000, 864.3494],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_hkl(ctax):
    hkl = (0, 0, 3)
    R0 = False
    rez = ctax.cooper_nathans(hkl=hkl, en=0, R0=R0)
    mat = np.array(
        [
            [33305.0843, 33224.4963, -2651.8290, -5152.9962],
            [33224.4963, 33305.2609, -2651.8526, -5153.0102],
            [-2651.8290, -2651.8526, 1983.2037, 448.8024],
            [-5152.9962, -5153.0102, 448.8024, 864.3494],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_projection(ctax):
    hkl = (0, 0, 3)
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = False
    rez = ctax.cooper_nathans(hkl=hkl, en=0, projection=projection, R0=R0)
    mat = np.array(
        [
            [1.3306e05, -5.3037e03, -1.7660e-01, -1.0306e04],
            [-5.3037e03, 1.9832e03, 2.3558e-02, 4.4880e02],
            [-1.7660e-01, 2.3558e-02, 1.6135e02, 1.4003e-02],
            [-1.0306e04, 4.4880e02, 1.4003e-02, 8.6435e02],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_list(ctax):
    R0 = False
    rez_list = ctax.cooper_nathans(hkl=[(0, 0, 3), (0, 0, -3)], en=[0, 1], projection=None, R0=R0)
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
def ctax():
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax = TAS(fixed_ef=4.8)
    ctax.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Sample.from_json(sample_json_path)
    ctax.mount_sample(sample)

    return ctax
