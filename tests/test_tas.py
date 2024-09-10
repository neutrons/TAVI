import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample.xtal import Xtal
from tavi.utilities import MotorAngles, Peak


def test_find_two_theta():
    ctax = TAS()
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Xtal.from_json("./test_data/test_samples/nitio3.json")
    ctax.mount_sample(nitio3)

    two_theta = ctax.calculate_two_theta(hkl=(0, 0, 3), ei=4.799999)
    assert np.allclose(two_theta, 53.240000, atol=1e-1)
    two_theta = ctax.calculate_two_theta(hkl=(0.5, 0.5, 0), ei=4.799999)
    assert np.allclose(two_theta, 48.489200, atol=1e-1)


def test_calc_ub_from_2_peaks_hb3():
    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.942840, 0.014534, -0.332928]

    tas = TAS()
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)
    tas.mount_sample(Xtal(lattice_params))

    ei = 13.500172
    ef = 13.505137
    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak(hkl=(0, 0, 2), angles=angles1, ei=ei, ef=ef)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak(hkl=(0, 2, 0), angles=angles2, ei=ei, ef=ef)

    tas.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(tas.sample.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(tas.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
    assert np.allclose(tas.sample.u, u, atol=1e-2)
    assert np.allclose(tas.sample.v, v, atol=1e-2)

    angles_1 = tas.find_angles(peak=(0, 0, 2), ei=13.500172)
    assert np.allclose(angles_1, angles1, atol=1e-2)

    # angles = tas.find_angles(peak=(0, 2, 0), ei=13.500172, ef=13.505137)
    # assert np.allclose(angles_list[1], angles, atol=1e-1)

    # swap peaks and calculate again
    tas.calculate_ub_matrix(peaks=(peak2, peak1))

    assert np.allclose(tas.sample.ub_mat, ub_matrix, atol=1e-2)
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
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Xtal.from_json("./test_data/test_samples/nitio3.json")

    ctax.mount_sample(nitio3)

    ei = 4.799999
    ef = 4.799998
    angles1 = MotorAngles(two_theta=53.240000, omega=32.865000, sgl=2.310799, sgu=2.021008)
    peak1 = Peak(hkl=(0, 0, 3), angles=angles1, ei=ei, ef=ef)
    angles2 = MotorAngles(two_theta=48.489200, omega=-59.507500, sgl=2.310799, sgu=2.021008)
    peak2 = Peak(hkl=(0.5, 0.5, 0), angles=angles2, ei=ei, ef=ef)

    ctax.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(ctax.sample.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(ctax.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(ctax.sample.in_plane_ref, in_plane_ref, atol=1e-2)

    # angles = ctax.find_angles(peak=(0, 0, 3), ei=4.799999, ef=4.799998)
    # assert np.allclose(angles_list[0], angles, atol=1e-1)
    # angles = ctax.find_angles(peak=(0.5, 0.5, 0), ei=4.799999, ef=4.799998)
    # assert np.allclose(angles_list[1], angles, atol=1e-1)


def test_calc_ub_from_2_peaks_hb1():
    hb1 = TAS()
    hb1_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    hb1.load_instrument_params_from_json(hb1_json)

    lattice_params = (3.939520, 3.939520, 3.941957, 90.0, 90.0, 90.0)
    xtal = Xtal(lattice_params=lattice_params)
    hb1.mount_sample(xtal)
    ub_matrix = np.array(
        [
            [0.181283, 0.177659, -0.002718],
            [0.177656, -0.181301, -0.001380],
            [-0.002909, -0.000917, -0.253663],
        ]
    )
    plane_normal = [0.010097, 0.999933, -0.005550]
    in_plane_ref = [0.999892, -0.010155, -0.010658]

    ei = 13.499993
    ef = 13.506112
    angles1 = MotorAngles(two_theta=-52.449830, omega=-26.825625, sgl=-0.578500, sgu=0.318000)
    peak1 = Peak(hkl=(1, 1, 0), angles=angles1, ei=ei, ef=ef)
    angles2 = MotorAngles(two_theta=-36.368334, omega=-108.827250, sgl=-0.578500, sgu=0.318000)
    peak2 = Peak(hkl=(0, 0, 1), angles=angles2, ei=ei, ef=ef)

    hb1.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(hb1.sample.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(hb1.sample.in_plane_ref, in_plane_ref, atol=1e-2)

    # angles = hb1.find_angles(peak=(1, 1, 0), ei=13.499993, ef=13.506112)
    # assert np.allclose(angles_list[0], angles, atol=1e-1)
    # angles = hb1.find_angles(peak=(0, 0, 1), ei=13.499993, ef=13.506112)
    # assert np.allclose(angles_list[1], angles, atol=1e-1)

    # swap peaks and calculate again
    hb1.calculate_ub_matrix(peaks=(peak2, peak1))

    assert np.allclose(hb1.sample.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.plane_normal, plane_normal, atol=1e-2)
    # assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
