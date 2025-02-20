import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def test_calc_ub_from_2_peaks_takin():

    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    spice_ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [-0.164330, -0.304247, 0.058788],
            [0.272815, -0.013290, 0.002566],
        ]
    )

    # u = [0.15623, 2.83819, -1.88465]
    # v = [-0.00060, 1.03219, 5.33915]

    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.942840, 0.014534, -0.332928]

    ei = 13.500172
    ef = 13.505137

    tas = TAS(fixed_ef=ef, fixed_ei=ei, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    tas.mount_sample(Sample(lattice_params))

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ub_conf = tas.calculate_ub_matrix(peaks=(peak1, peak2))
    # u_cal, v_cal = tas.sample.ub_matrix_to_uv(tas.sample.ub_mat)

    assert np.allclose(ub_conf.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(ub_conf.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(ub_conf.in_plane_ref, in_plane_ref, atol=1e-2)
    # assert np.allclose(u_cal, u, atol=1e-2)
    # assert np.allclose(v_cal, v, atol=1e-2)

    # u_cal, v_cal = tas.sample.ub_matrix_to_uv(spice_to_mantid(spice_ub_matrix))
    # assert np.allclose(u_cal, u, atol=1e-2)
    # assert np.allclose(v_cal, v, atol=1e-2)

    angles_1 = tas.calculate_motor_angles((0, 0, 2))
    for name, value in angles_1._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-1)

    angles_2 = tas.calculate_motor_angles((0, 2, 0))
    for name, value in angles_2._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e0)

    # swap peaks and calculate again
    ub_conf = tas.calculate_ub_matrix(peaks=(peak2, peak1))
    #  u_cal, v_cal = tas.sample.ub_matrix_to_uv(tas.sample.ub_mat)

    assert np.allclose(ub_conf.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(ub_conf.plane_normal, plane_normal, atol=1e-2)
    # assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
    # assert np.allclose(u_cal, u, atol=1e-2)
    # assert np.allclose(v_cal, v, atol=1e-2)


def test_calc_ub_from_2_peaks_hb3():
    lattice_params = (4.128474, 4.128474, 6.651507, 90, 90, 120)
    ub_matrix = np.array(
        [
            [-0.000328, -0.004396, -0.150319],
            [-0.279308, -0.152332, 0.000314],
            [-0.014647, 0.234528, -0.002614],
        ]
    )

    plane_normal = [-0.017470, -0.052341, 0.998476]
    in_plane_ref = [-0.999847, -0.002086, -0.017385]

    tas = TAS(fixed_ef=14.7, spice_convention=True)
    hb3_json = "./src/tavi/instrument/instrument_params/hb3_mnte.json"
    tas.load_instrument_params_from_json(hb3_json)
    tas.mount_sample(Sample(lattice_params))

    angles1 = MotorAngles(two_theta=41.545383, omega=20.840750, sgl=1.001000, sgu=-3.000750)
    peak1 = Peak(hkl=(0, 0, 2), angles=angles1)
    angles2 = MotorAngles(two_theta=38.526091, omega=-70.614125, sgl=1.001000, sgu=-3.000750)
    peak2 = Peak(hkl=(1, 0, 0), angles=angles2)

    ub_conf = tas.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(ub_conf.ub_mat, ub_matrix, atol=1e-1)
    assert np.allclose(ub_conf.plane_normal, plane_normal, atol=1e-1)
    assert np.allclose(ub_conf.in_plane_ref, in_plane_ref, atol=1e-1)

    angles_1 = tas.calculate_motor_angles(hkl=(0, 0, 2))
    for name, value in angles_1._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-1)

    angles_2 = tas.calculate_motor_angles(hkl=(1, 0, 0))
    for name, value in angles_2._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e-1)

    assert np.allclose(tas.sample.ub_conf.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(tas.sample.ub_conf.in_plane_ref, in_plane_ref, atol=1e-2)


def test_calc_ub_from_2_peaks_ctax():
    ub_matrix = np.array(
        [
            [-0.016934, -0.026164, -0.071871],
            [-0.20102, -0.192954, 0.007764],
            [-0.108217, 0.12038, -0.003176],
        ]
    )
    plane_normal = [-0.04032, 0.035237, 0.998565]
    in_plane_ref = [-0.993257, 0.107299, -0.043892]

    ctax = TAS(fixed_ef=4.799998, spice_convention=True)
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Sample.from_json("./test_data/test_samples/nitio3.json")

    ctax.mount_sample(nitio3)

    angles1 = MotorAngles(two_theta=53.240000, omega=32.865000, sgl=2.310799, sgu=2.021008)
    peak1 = Peak(hkl=(0, 0, 3), angles=angles1)
    angles2 = MotorAngles(two_theta=48.489200, omega=-59.507500, sgl=2.310799, sgu=2.021008)
    peak2 = Peak(hkl=(0.5, 0.5, 0), angles=angles2)

    ctax.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(ctax.sample.ub_conf.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(ctax.sample.ub_conf.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(ctax.sample.ub_conf.in_plane_ref, in_plane_ref, atol=1e-2)

    angles = ctax.calculate_motor_angles(hkl=(0, 0, 3))
    for name, value in angles._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-1)
    angles = ctax.calculate_motor_angles(hkl=(0.5, 0.5, 0))
    for name, value in angles._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e-1)


def test_calc_ub_from_2_peaks_hb1():
    ei = 13.499993
    ef = 13.506112
    hb1 = TAS(fixed_ef=ef, fixed_ei=ei, spice_convention=False)
    hb1_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    hb1.load_instrument_params_from_json(hb1_json)

    lattice_params = (3.939520, 3.939520, 3.941957, 90.0, 90.0, 90.0)
    xtal = Sample(lattice_params=lattice_params)
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

    angles1 = MotorAngles(two_theta=-52.449830, omega=-26.825625, sgl=-0.578500, sgu=0.318000)
    peak1 = Peak(hkl=(1, 1, 0), angles=angles1)
    angles2 = MotorAngles(two_theta=-36.368334, omega=-108.827250, sgl=-0.578500, sgu=0.318000)
    peak2 = Peak(hkl=(0, 0, 1), angles=angles2)

    hb1.calculate_ub_matrix(peaks=(peak1, peak2))

    assert np.allclose(hb1.sample.ub_conf.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.ub_conf.plane_normal, plane_normal, atol=1e-2)
    assert np.allclose(hb1.sample.ub_conf.in_plane_ref, in_plane_ref, atol=1e-2)

    angles = hb1.calculate_motor_angles(hkl=(1, 1, 0))
    for name, value in angles._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-1)
    angles = hb1.calculate_motor_angles(hkl=(0, 0, 1))
    for name, value in angles._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e-1)

    # swap peaks and calculate again
    hb1.calculate_ub_matrix(peaks=(peak2, peak1))

    assert np.allclose(hb1.sample.ub_conf.ub_mat, ub_matrix, atol=1e-2)
    assert np.allclose(hb1.sample.ub_conf.plane_normal, plane_normal, atol=1e-2)
    # assert np.allclose(tas.sample.in_plane_ref, in_plane_ref, atol=1e-2)
