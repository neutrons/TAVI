import json

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import Normalize
from mpl_toolkits.axisartist import Axes

from tavi.data.scan import Scan
from tavi.data.scan_group import ScanGroup
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

# plt.rcParams.update({"font.size": 8})


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

    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")
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
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")
    with pytest.raises(AttributeError) as e_info:
        rez = tas.cooper_nathans(hkle=hkle, axes=projection)
    assert "sample is missing in TAS(fixed_ei=5, fixed_ef=None, convention=Spice)." in str(e_info.value)


def test_out_of_reach(sample):
    tas = TAS(fixed_ei=5)
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (10, 10, 10, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")

    rez = tas.cooper_nathans(hkle=hkle, axes=projection)
    assert rez is None
    assert "Cannot get two_theta for hkl=(10, 10, 10), ei=5 meV, ef=5 meV. Triangle cannot be closed." in tas.err_msg


def test_reso_str(sample):
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"

    tas = TAS(fixed_ei=5)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (0.1, 0.1, 0, 0)
    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")
    rez = tas.cooper_nathans(hkle=hkle, axes=projection)

    assert len(str(rez).split("\n")) == 31


def test_plot_ellipses(sample):
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"

    tas = TAS(fixed_ei=5)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample)

    hkle = (0.1, 0.1, 0, 0)
    projection = ((1, 1, 0), (0, 1, 0), (0, 0, 1), "en")
    rez = tas.cooper_nathans(hkle=hkle, axes=projection)

    # plotting
    fig = plt.figure(figsize=(10, 6), constrained_layout=True)
    rez.plot_ellipses(fig)
    plt.show()


def test_generate_hkle_list():
    axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
    hkle = TAS.generate_hkle(
        axes=axes,
        grid=((-0.5, 0.15, 0.05), 3, 0, (0, 4.1, 0.4)),
    )
    assert len(hkle) == 143


def test_plot_ellipsoids_contour():
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"

    param_dict = {"fixed_ef": 4.8}
    tas = TAS(**param_dict)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(Sample.from_json("./test_data/test_samples/nitio3.json"))

    # calculate resolution ellipses
    axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
    grid = ((-0.5, 5.5, 0.05), 3, 0, (0, 4.1, 0.4))

    axes = ("en", (1, 1, 0), (0, 0, 1), (1, -1, 0))
    grid = ((0, 4.1, 0.4), (-0.5, 5.5, 0.05), 3, 0)

    hkle_list = tas.generate_hkle(grid, axes=axes)
    rez_list = tas.cooper_nathans(hkle=hkle_list, axes=axes)

    assert len(rez_list) == 1320

    # generate plot
    p = Plot2D()
    axes = [i for i, val in enumerate(grid) if isinstance(val, tuple) and len(val) == 3]
    for rez in filter(None, rez_list):
        e_inco = rez.get_ellipse(axes=axes, PROJECTION=True)
        e_co = rez.get_ellipse(axes=axes, PROJECTION=False)
        p.add_reso(e_co, c="k", linestyle="solid")
        p.add_reso(e_inco, c="k", linestyle="dashed")

    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)
    p.plot(ax)
    plt.show()


def test_plot_data_contour():
    path_to_spice_folder = "test_data/exp424/"
    scan_nums = "42-48, 70-75"

    scan_data_2d_params = {
        "axes": ("en", "qh", "detector"),
        "norm_to": (1, "mcu"),
        "grid": ((0, 4.1, 0.1), 0.025),
    }

    contour_params = {
        "cmap": "turbo",
        "norm": Normalize(vmin=0, vmax=1),
        # "norm": LogNorm(vmin=1e-1, vmax=1e4),
    }

    scan_list = []
    for part in scan_nums.split(","):
        start, end = map(int, part.strip().split("-"))
        scan_list.extend(range(start, end + 1))
    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for num in scan_list]
    sg = ScanGroup(scans)

    scan_data_2d = sg.combine_data(**scan_data_2d_params)
    p = Plot2D()

    p.add_contour(scan_data_2d, **contour_params)
    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)
    im = p.plot(ax)
    fig.colorbar(im, ax=ax)

    plt.show()


def test_plot_ellipsoids_contour_oplot():
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"

    param_dict = {"fixed_ef": 4.8}
    tas = TAS(**param_dict)
    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(Sample.from_json("./test_data/test_samples/nitio3.json"))

    # calculate resolution ellipses
    axes = ("en", (1, 1, 0), (0, 0, 1), (1, -1, 0))
    grid = ((0, 4.1, 0.4), (-0.5, 0.15, 0.05), 3, 0)

    hkle_list = tas.generate_hkle(grid, axes)
    rez_list = tas.cooper_nathans(hkle=hkle_list, axes=axes)

    # generate data
    path_to_spice_folder = "test_data/exp424/"
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for num in scan_list]
    sg = ScanGroup(scans)

    scan_data_2d_params = {
        "axes": ("en", "qh", "detector"),
        "norm_to": (1, "mcu"),
        "grid": ((0, 4.5, 0.1), 0.025),
    }

    scan_data_2d = sg.combine_data(**scan_data_2d_params)

    # generate plot
    p = Plot2D()

    contour_params = {
        "cmap": "turbo",
        "vmax": 1,
        "vmin": 0,
    }

    p.add_contour(scan_data_2d, **contour_params)

    axes_ellip = [i for i, val in enumerate(grid) if isinstance(val, tuple) and len(val) == 3]
    for rez in filter(None, rez_list):
        e_co = rez.get_ellipse(axes=axes_ellip, PROJECTION=False)
        e_inco = rez.get_ellipse(axes=axes_ellip, PROJECTION=True)
        p.add_reso(e_co, c="k", linestyle="solid")
        p.add_reso(e_inco, c="k", linestyle="dashed")

    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)
    im = p.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()
