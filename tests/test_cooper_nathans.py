import matplotlib.pyplot as plt
import numpy as np
import pytest

from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import (
    UBConf,
    b_mat_from_ub_matrix,
    mantid_to_spice,
    plane_normal_from_two_peaks,
    u_mat_from_ub_matrix,
    uv_to_ub_matrix,
)


def test_validate_instrument_parameter(ctax):
    cn = CooperNathans(instrument=ctax)
    cn.validate_instrument_parameters()


def test_generate_hkle_list(ctax):
    cn = CooperNathans(instrument=ctax)
    hkleief_list = cn.generate_hkleief_list(
        hkle=[(h, h, 3, en) for h in np.arange(0, 1, 0.1) for en in np.arange(0, 5, 0.5)],
    )
    assert len(hkleief_list) == 100
    (h, k, l), ei, ef = hkleief_list[0]
    assert np.allclose((h, k, l), (0, 0, 3))
    assert np.allclose((ei, ef), (4.8, 4.8))


def test_local_q(ctax):
    hkle = (0, 0, 3, 0)
    rez = ctax.cooper_nathans(hkle=hkle, axes=None)
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
    hkle = (0, 0, 3, 0)
    rez = ctax.cooper_nathans(hkle=hkle)
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
    hkle = (0, 0, 3, 0)
    axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
    rez = ctax.cooper_nathans(hkle=hkle, axes=axes)
    mat = np.array(
        [
            [1.3306e05, -5.3037e03, -1.7660e-01, -1.0306e04],
            [-5.3037e03, 1.9832e03, 2.3558e-02, 4.4880e02],
            [-1.7660e-01, 2.3558e-02, 1.6135e02, 1.4003e-02],
            [-1.0306e04, 4.4880e02, 1.4003e-02, 8.6435e02],
        ]
    )
    assert np.allclose(rez.mat, mat, atol=1e-1)


def test_projection_out_of_reach():
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = TAS(fixed_ei=4.8)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(Sample.from_json("./test_data/test_samples/nitio3.json"))

    hkle = (-1 / 2, -1 / 2, 3, 3.2)
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    rez = tas.cooper_nathans(hkle=hkle, axes=projection)
    assert not rez
    assert "Cannot get theta_a for ef=1.6 meV. Bragg condition cannot be fulfilled." in tas.err_msg


def test_list(ctax):
    rez_list = ctax.cooper_nathans(
        hkle=[
            (0, 0, 3, 0),
            (0, 0, -3, 0),
            (0, 0, 3, 1),
            (0, 0, -3, 1),
        ],
        axes=None,
    )
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


@pytest.fixture
def rescal_params():
    instrument_params = {
        "monochromator": {
            "type": "PG002",
            "mosaic_h": 30,
            "mosaic_v": 30,
            "sense": "+",
        },
        "goniometer": {"sense": "-", "type": "Y,-Z,X"},
        "analyzer": {
            "type": "Pg002",
            "mosaic_h": 30,
            "mosaic_v": 30,
            "sense": "+",
        },
        "collimators": {
            "h_pre_mono": 120,
            "h_pre_sample": 60,
            "h_post_sample": 60,
            "h_post_ana": 60,
            "v_pre_mono": 120,
            "v_pre_sample": 120,
            "v_post_sample": 120,
            "v_post_ana": 120,
        },
    }

    lattice_params = (6.28319, 6.28319, 6.28319, 90, 90, 90)
    sample = Sample(lattice_params)
    sample.set_mosaic(10, 10)
    u = (1, 0, 0)
    v = (0, 1, 0)
    ub_matrix_mantid = uv_to_ub_matrix(u, v, lattice_params)
    plane_normal_mantid, in_plane_ref_mantid = plane_normal_from_two_peaks(
        u_mat_from_ub_matrix(ub_matrix_mantid),
        b_mat_from_ub_matrix(ub_matrix_mantid),
        (1, 0, 0),
        (0, 1, 0),
    )

    sample.ub_conf = UBConf(
        ub_mat=mantid_to_spice(ub_matrix_mantid),
        plane_normal=mantid_to_spice(plane_normal_mantid),
        in_plane_ref=mantid_to_spice(in_plane_ref_mantid),
    )
    return instrument_params, sample


def test_rescal_comparison(rescal_params):
    instrument_params, sample = rescal_params
    tas = TAS(fixed_ei=14.7)
    tas._load_instrument_parameters(instrument_params)
    tas.mount_sample(sample)

    rez_hkl = tas.cooper_nathans(hkle=(1, 2, 0, 0))
    rez_hkl.plot()

    rez_invA = tas.cooper_nathans(hkle=(1, 2, 0, 0), axes=None)
    assert np.allclose(
        rez_invA.mat,
        [
            [6654, -1555, 0, -179.7],
            [-1555, 3.104e04, 0, 2963],
            [0, 0, 628.9, 0],
            [-179.7, 2963, 0, 285.1],
        ],
        rtol=0.01,
    )

    plt.show()


@pytest.fixture
def takin_params():
    instrument_params = {
        "monochromator": {
            "type": "PG002",
            "mosaic_h": 30,
            "mosaic_v": 30,
            "sense": "-",
        },
        "goniometer": {"sense": "+", "type": "Y,-Z,X"},
        "analyzer": {
            "type": "Pg002",
            "mosaic_h": 30,
            "mosaic_v": 30,
            "sense": "-",
        },
        "collimators": {
            "h_pre_mono": 40,
            "h_pre_sample": 40,
            "h_post_sample": 40,
            "h_post_ana": 80,
            "v_pre_mono": 600,
            "v_pre_sample": 600,
            "v_post_sample": 600,
            "v_post_ana": 600,
        },
    }

    lattice_params = (5, 5, 5, 90, 90, 90)
    sample = Sample(lattice_params)
    sample.set_mosaic(0, 0)
    u = (1.067 / 2, 0.9995, 0)
    v = (-1.067 / 2, 0.9995, 0)
    ub_matrix_mantid = uv_to_ub_matrix(u, v, lattice_params)
    plane_normal_mantid, in_plane_ref_mantid = plane_normal_from_two_peaks(
        u_mat_from_ub_matrix(ub_matrix_mantid),
        b_mat_from_ub_matrix(ub_matrix_mantid),
        (1, 0, 0),
        (0, 1, 0),
    )

    sample.ub_conf = UBConf(
        ub_mat=mantid_to_spice(ub_matrix_mantid),
        plane_normal=mantid_to_spice(plane_normal_mantid),
        in_plane_ref=mantid_to_spice(in_plane_ref_mantid),
    )
    return instrument_params, sample


def test_takin_comparison(takin_params):
    instrument_params, sample = takin_params
    tas = TAS(fixed_ef=2.661513)
    tas._load_instrument_parameters(instrument_params)
    tas.mount_sample(sample)

    rez = tas.cooper_nathans(hkle=(0.8488, 0, 0, 0), axes=None)

    assert np.allclose(rez.coh_fwhms(0), 0.0112653, rtol=0.01)
    assert np.allclose(rez.coh_fwhms(1), 0.00486532, rtol=0.01)
    assert np.allclose(rez.coh_fwhms(2), 0.19814, rtol=0.01)
    assert np.allclose(rez.coh_fwhms(3), 0.0149406, rtol=0.01)

    assert np.allclose(rez.incoh_fwhms(0), 0.0112775, rtol=0.01)
    assert np.allclose(rez.incoh_fwhms(1), 0.0164951, rtol=0.01)
    assert np.allclose(rez.incoh_fwhms(2), 0.19814, rtol=0.01)
    assert np.allclose(rez.incoh_fwhms(3), 0.0506693, rtol=0.01)

    assert np.allclose(
        rez.mat,
        [
            [43694.7, 2303.5, 0, -1109.74],
            [2303.5, 234257, 0, -72887.1],
            [0, 0, 141.244, 0],
            [-1109.74, -72887.1, 0, 24841.6],
        ],
        rtol=0.01,
    )

    rez = tas.cooper_nathans(hkle=(0.8488, 0, 0, 0))
    assert np.allclose(
        rez.mat,
        [
            [68999.9, 3637.54, 0, -1394.54],
            [3637.54, 369923, 0, -91592.6],
            [0, 0, 223.044, 0],
            [-1394.54, -91592.6, 0, 24841.6],
        ],
        rtol=0.01,
    )
    rez.plot()
    plt.show()
