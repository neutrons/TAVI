import json

import matplotlib.pyplot as plt
import numpy as np
import pytest
from mpl_toolkits.axisartist import Axes

from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample
from tavi.ub_algorithm import (
    UBConf,
    b_mat_from_ub_matrix,
    mantid_to_spice,
    plane_normal_from_two_peaks,
    spice_to_mantid,
    u_mat_from_ub_matrix,
    uv_to_ub_matrix,
)


@pytest.fixture
def sample():
    lattice_params = (10, 10, 10, 90, 90, 90)
    sample = Sample(lattice_params)
    sample.set_mosaic(30, 30)
    u = (1, 0, 0)
    v = (0, 1, 0)
    ub_matrix_mantid = uv_to_ub_matrix(u, v, lattice_params)
    ub_matrix_spice = mantid_to_spice(ub_matrix_mantid)

    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    b_mat = b_mat_from_ub_matrix(ub_matrix_mantid)
    u_mat_mantid = u_mat_from_ub_matrix(ub_matrix_mantid)
    plane_normal_mantid, in_plane_ref_mantid = plane_normal_from_two_peaks(
        u_mat_mantid, b_mat, projection[0], projection[1]
    )

    sample.ub_conf = UBConf(
        ub_mat=ub_matrix_spice,
        plane_normal=mantid_to_spice(plane_normal_mantid),
        in_plane_ref=mantid_to_spice(in_plane_ref_mantid),
    )
    return sample


def test_uv_to_ub_matrix():
    u = (1, 0, 0)
    v = (0, 1, 0)
    lattice_params = (10, 10, 10, 90, 90, 90)

    ub_matrix_calc = uv_to_ub_matrix(u, v, lattice_params)
    ub_matrix_calc = mantid_to_spice(ub_matrix_calc)
    ub_spice = ((0, 0.1, 0), (-0.1, 0, 0), (0, 0, 0.1))

    assert np.allclose(ub_matrix_calc, ub_spice, atol=1e-2)


def test_b_mat_from_ub_matrix():
    lattice_params = (10, 10, 10, 90, 90, 90)
    sample = Sample(lattice_params)

    ub_spice = ((0, 0.1, 0), (-0.1, 0, 0), (0, 0, 0.1))
    b_mat_ub = b_mat_from_ub_matrix(spice_to_mantid(ub_spice))

    assert np.allclose(b_mat_ub, sample.b_mat)


def test_u_mat_from_ub_matrix():
    ub_spice = ((0, 0.1, 0), (-0.1, 0, 0), (0, 0, 0.1))
    u_mat_cal = u_mat_from_ub_matrix(spice_to_mantid(ub_spice))
    u_mat = ((0, 1, 0), (0, 0, 1), (1, 0, 0))
    assert np.allclose(u_mat_cal, u_mat)


def test_uv_to_plane_normal():
    ub_spice = ((0, 0.1, 0), (-0.1, 0, 0), (0, 0, 0.1))

    b_mat = b_mat_from_ub_matrix(spice_to_mantid(ub_spice))
    u_mat = u_mat_from_ub_matrix(spice_to_mantid(ub_spice))

    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    plane_normal, in_plane_ref = plane_normal_from_two_peaks(u_mat, b_mat, projection[0], projection[1])
    assert np.allclose(plane_normal, (0, 1, 0))
    assert np.allclose(in_plane_ref, (0, 0, 1))


def test_instrument_sample_setup(sample):
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"
    with open(instrument_config_json_path, "r", encoding="utf-8") as file:
        config_params_dict = json.load(file)

    tas = TAS(fixed_ei=5)
    tas._load_instrument_parameters(config_params_dict)
    assert np.allclose(tas.fixed_ei, 5.0)

    tas.mount_sample(sample)

    lattice_params = (10, 10, 10, 90, 90, 90)
    assert np.allclose(tas.sample.lattice_params, lattice_params)


def test_instrument_error_handling():
    tas = TAS(fixed_ei=5)
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"
    tas.load_instrument_params_from_json(instrument_config_json_path)

    hkle = (0.1, 0.1, 0, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    with pytest.raises(AttributeError) as e_info:
        rez = tas.cooper_nathans(hkle=hkle, projection=projection)
    assert "sample is missing in TAS(fixed_ei=5, fixed_ef=None, convention=Spice)." in str(e_info.value)


def test_out_of_reach(sample):
    tas = TAS(fixed_ei=5)
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (10, 10, 10, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

    rez = tas.cooper_nathans(hkle=hkle, projection=projection)
    assert rez is None
    assert (
        "Cannot get two_theta for hkl=(10, 10, 10), ei=5.00 meV, ef=5.00 meV. Triangle cannot be closed." in tas.err_msg
    )


def test_reso_str(sample):
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"

    tas = TAS(fixed_ei=5)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (0.1, 0.1, 0, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    # projection = None
    rez = tas.cooper_nathans(hkle=hkle, projection=projection)

    assert len(str(rez).split("\n")) == 31


def test_plot_ellipses(sample):
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"

    tas = TAS(fixed_ei=5)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (0.1, 0.1, 0, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    # projection = None
    rez = tas.cooper_nathans(hkle=hkle, projection=projection)

    # plotting
    fig = plt.figure(figsize=(10, 6))

    for i, indices in enumerate([(0, 3), (1, 3), (2, 3), (0, 1), (1, 2), (0, 2)]):
        ellipse_co = rez.get_ellipse(axes=indices, PROJECTION=False)
        ellipse_inco = rez.get_ellipse(axes=indices, PROJECTION=True)

        p = Plot2D()
        if indices == (2, 3):
            p.add_reso(ellipse_co, c="k", linestyle="solid", label="Coherent")
            p.add_reso(ellipse_inco, c="k", linestyle="dashed", label="Incoherent")

        else:
            p.add_reso(ellipse_co, c="k", linestyle="solid")
            p.add_reso(ellipse_inco, c="k", linestyle="dashed")

        ax = fig.add_subplot(
            int(f"23{i + 1}"),
            axes_class=Axes,
            grid_helper=p.grid_helper(ellipse_co.angle),
        )
        p.plot(ax)

    fig.tight_layout(pad=2)
    plt.show()
