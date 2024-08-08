import numpy as np
from tavi.instrument.tas import TAS
from tavi.sample.xtal import Xtal


def test_load_from_json():
    tax = TAS()
    cg4c = "./src/tavi/instrument/instrument_params/cg4c.json"
    nitio3 = "./test_data/test_samples/nitio3.json"
    tax.load_instrument_from_json(cg4c)
    tax.load_sample_from_json(nitio3)

    assert tax.analyzer.type == "Pg002"
    assert np.allclose(tax.analyzer.d_spacing, 3.35416)
    assert np.allclose(tax.sample.a, 5.034785)
    assert np.allclose(tax.sample.gamma, 120)


def test_calc_ub_from_2_peaks_hb3():
    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    b_matrix = np.array(
        [
            [0.3230, 0.1615, 0.0000],
            [0.0000, 0.2797, -0.0000],
            [0.0000, 0.0000, 0.1766],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.942840, 0.014534, -0.332928]

    tas = TAS()
    takin = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_from_json(takin)
    tas.load_sample(Xtal(lattice_params))

    peak_list = [
        (0, 0, 2),
        (0, 2, 0),
    ]
    # (s2, s1, sgl, sgu)
    angles_list = [
        (-51.530388, -45.220125, -0.000500, -2.501000),
        (-105.358735, 17.790125, -0.000500, -2.501000),
    ]

    tas.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)

    assert np.allclose(tas.sample.ub_matrix, ub_matrix, atol=1e-2)
    assert np.allclose(tas.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
    assert np.allclose(tas.sample.u, u, atol=1e-2)
    assert np.allclose(tas.sample.v, v, atol=1e-2)

    angles = tas.find_angles(peak=(0, 0, 2), ei=13.500172, ef=13.505137)
    assert np.allclose(angles_list[0], angles, atol=1e-2)

    angles = tas.find_angles(peak=(0, 2, 0), ei=13.500172, ef=13.505137)
    assert np.allclose(angles_list[1], angles, atol=1e-1)

    # swap peaks and calculate again
    tas.find_ub(
        peaks=(peak_list[1], peak_list[0]),
        angles=(angles_list[1], angles_list[0]),
        ei=13.500172,
        ef=13.505137,
    )
    assert np.allclose(tas.sample.ub_matrix, ub_matrix, atol=1e-2)
    assert np.allclose(tas.sample.plane_normal, plane_normal, atol=1e-2)
    # assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
    assert np.allclose(tas.sample.u, u, atol=1e-2)
    assert np.allclose(tas.sample.v, v, atol=1e-2)


def test_calc_ub_from_2_peaks_ctax():
    ub_matrix = np.array(
        [
            [-0.016934, -0.026164, -0.071871],
            [-0.108217, 0.12038, -0.003176],
            [0.20102, 0.192954, -0.007764],
        ]
    )
    plane_normal = [-0.04032, 0.998565, -0.035237]
    in_plane_ref = [-0.993257, -0.043892, -0.107299]

    ctax = TAS()
    cg4c = "./src/tavi/instrument/instrument_params/cg4c.json"
    nitio3 = "./test_data/test_samples/nitio3.json"
    ctax.load_instrument_from_json(cg4c)
    ctax.load_sample_from_json(nitio3)

    peak_list = [
        (0, 0, 3),
        (0.5, 0.5, 0),
    ]
    # (s2, s1, sgl, sgu)
    angles_list = [
        (53.240000, 32.865000, 2.310799, 2.021008),
        (48.489200, -59.507500, 2.310799, 2.021008),
    ]

    ctax.find_ub(peaks=peak_list, angles=angles_list, ei=4.799999, ef=4.799998)

    assert np.allclose(ctax.sample.ub_matrix, ub_matrix, atol=1e-2)
    assert np.allclose(ctax.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(ctax.sample.in_plane_ref, in_plane_ref, atol=1e-2)

    angles = ctax.find_angles(peak=(0, 0, 3), ei=4.799999, ef=4.799998)
    assert np.allclose(angles_list[0], angles, atol=1e-1)
    angles = ctax.find_angles(peak=(0.5, 0.5, 0), ei=4.799999, ef=4.799998)
    assert np.allclose(angles_list[1], angles, atol=1e-1)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)


def test_calc_ub_from_2_peaks_hb1():
    hb1 = TAS()
    hb1.load_instrument_from_json("./src/tavi/instrument/instrument_params/takin_test.json")

    lattice_params = (3.939520, 3.939520, 3.941957, 90.0, 90.0, 90.0)

    hb1.load_sample(Xtal(lattice_params=lattice_params))
    ub_matrix = np.array(
        [
            [0.181283, 0.177659, -0.002718],
            [0.177656, -0.181301, -0.001380],
            [-0.002909, -0.000917, -0.253663],
        ]
    )
    plane_normal = [0.010097, 0.999933, -0.005550]
    in_plane_ref = [0.999892, -0.010155, -0.010658]

    peak_list = [(1, 1, 0), (0, 0, 1)]
    # (s2, s1, sgl, sgu)
    angles_list = [
        (-52.449830, -26.825625, -0.578500, 0.318000),
        (-36.368334, -108.827250, -0.578500, 0.318000),
    ]

    hb1.find_ub(peaks=peak_list, angles=angles_list, ei=13.499993, ef=13.506112)

    assert np.allclose(hb1.sample.ub_matrix, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(hb1.sample.in_plane_ref, in_plane_ref, atol=1e-2)

    angles = hb1.find_angles(peak=(1, 1, 0), ei=13.499993, ef=13.506112)
    assert np.allclose(angles_list[0], angles, atol=1e-1)
    angles = hb1.find_angles(peak=(0, 0, 1), ei=13.499993, ef=13.506112)
    assert np.allclose(angles_list[1], angles, atol=1e-1)

    # swap peaks and calculate again
    hb1.find_ub(
        peaks=(peak_list[1], peak_list[0]),
        angles=(angles_list[1], angles_list[0]),
        ei=13.499993,
        ef=13.506112,
    )

    assert np.allclose(hb1.sample.ub_matrix, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.plane_normal, plane_normal, atol=1e-2)
    # assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
