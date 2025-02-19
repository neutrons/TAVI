import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import r_matrix_with_minimal_tilt
from tavi.utilities import MotorAngles, Peak, UBConf


def test_find_two_theta():
    ctax = TAS(fixed_ef=4.799999)
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Sample.from_json("./test_data/test_samples/nitio3.json")
    ctax.mount_sample(nitio3)

    two_theta = ctax.get_two_theta(hkl=(0, 0, 3))
    assert np.allclose(two_theta, 53.240000, atol=1e-1)
    two_theta = ctax.get_two_theta(hkl=(0.5, 0.5, 0))
    assert np.allclose(two_theta, 48.489200, atol=1e-1)


def test_calculate_ub_matrix():
    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    sample = Sample(lattice_params)
    tas.mount_sample(sample)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ubconf = tas.calculate_ub_matrix((peak1, peak2))
    assert np.allclose(ubconf.ub_mat, ub_matrix, atol=1e-4)


def test_angles_to_r_mat():
    r_mat = np.array(
        [
            [7.04384932e-01, 3.09680706e-02, -7.09142332e-01],
            [8.72664626e-06, 9.99047460e-01, 4.36368240e-02],
            [7.09818194e-01, -3.07433097e-02, 7.03713707e-01],
        ]
    )
    angles = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    r_mat_cal = tas.goniometer.r_mat(angles)
    assert np.allclose(r_mat_cal, r_mat)


def test_r_mat_to_angles():
    r_mat = np.array(
        [
            [7.04384932e-01, 3.09680706e-02, -7.09142332e-01],
            [8.72664626e-06, 9.99047460e-01, 4.36368240e-02],
            [7.09818194e-01, -3.07433097e-02, 7.03713707e-01],
        ]
    )
    angles = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    angles_cal = tas.goniometer.angles_from_r_mat(r_mat, two_theta=angles.two_theta)
    for name, value in angles_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles, name), atol=1e-2)


def test_calculate_ub():

    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]
    sample = Sample(lattice_params)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ub_conf = tas.calculate_ub_matrix((peak1, peak2))
    assert np.allclose(ub_conf.ub_mat, ub_matrix, atol=1e-4)
    assert np.allclose(ub_conf.plane_normal, plane_normal, atol=1e-4)
    assert np.allclose(ub_conf.in_plane_ref, in_plane_ref, atol=1e-4)


def test_r_mat_with_minimal_tilt():

    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    sample = Sample(lattice_params)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    peak1 = Peak((0, 0, 2), angles1)
    angles2 = MotorAngles(two_theta=-105.358735, omega=17.790125, sgl=-0.000500, sgu=-2.501000)
    peak2 = Peak((0, 2, 0), angles2)

    ub_conf = tas.calculate_ub_matrix((peak1, peak2))

    r_mat = tas.goniometer.r_mat(angles1)
    r_mat_cal = r_matrix_with_minimal_tilt(
        hkl=(0, 0, 2), ei=13.505137, ef=13.505137, two_theta=angles1.two_theta, ub_conf=ub_conf
    )
    assert np.allclose(r_mat, r_mat_cal, atol=1e-3)


def test_calculate_motor_agnles():

    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    takin_json = "./src/tavi/instrument/instrument_params/takin_test.json"
    tas.load_instrument_params_from_json(takin_json)

    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.94290377, 0.01452569, -0.33274837]
    sample = Sample(lattice_params)
    sample.ub_conf = UBConf(ub_mat=ub_matrix, plane_normal=plane_normal, in_plane_ref=in_plane_ref)
    tas.mount_sample(sample)

    angels1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2), en=-0.005)
    angles1 = MotorAngles(two_theta=-51.530388, omega=-45.220125, sgl=-0.000500, sgu=-2.501000)
    for name, value in angels1_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-2)


def test_calculate_motor_angles_hb3():
    tas = TAS(fixed_ef=14.7)
    hb3_json = "./src/tavi/instrument/instrument_params/hb3_mnte.json"
    tas.load_instrument_params_from_json(hb3_json)

    lattice_params = (4.128474, 4.128474, 6.651507, 90, 90, 120)
    sample = Sample(lattice_params)
    ub_mat = np.array(
        [
            [-0.000328, -0.004396, -0.150319],
            [-0.279308, -0.152332, 0.000314],
            [-0.014647, 0.234528, -0.002614],
        ]
    )

    plane_normal = [-0.017470, -0.052341, 0.998476]
    in_plane_ref = [-0.999847, 0.002086, -0.017385]
    sample.ub_conf = UBConf(ub_mat=ub_mat, plane_normal=plane_normal, in_plane_ref=in_plane_ref)
    tas.mount_sample(sample)

    angles1 = MotorAngles(two_theta=41.545383, omega=20.840750, sgl=1.001000, sgu=-3.000750)
    r_mat = tas.goniometer.r_mat(angles1)
    r_mat_cal = r_matrix_with_minimal_tilt(
        hkl=(0, 0, 2), ei=14.7, ef=14.7, two_theta=angles1.two_theta, ub_conf=sample.ub_conf
    )
    assert np.allclose(r_mat, r_mat_cal, atol=1e-3)

    angels1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2))
    for name, value in angels1_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-2)

    angels2_cal = tas.calculate_motor_angles(hkl=(0, 0, 1))
    angles2 = MotorAngles(two_theta=38.526091, omega=-70.614125, sgl=1.001000, sgu=-3.000750)
    for name, value in angels2_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e-2)

    angels3_cal = tas.calculate_motor_angles(hkl=(1.3, 0, 1.3))
    angles3 = MotorAngles(two_theta=21.02, omega=6.71, sgl=1.001000, sgu=-3.000750)
    for name, value in angels3_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles3, name), atol=1e-2)
