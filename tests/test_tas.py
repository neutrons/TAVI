import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import r_matrix_with_minimal_tilt
from tavi.utilities import MotorAngles, Peak


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

    ub_cal = tas.calculate_ub_matrix((peak1, peak2))
    assert np.allclose(ub_cal, ub_matrix, atol=1e-4)


def test_r_mat():

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

    ub_cal = tas.calculate_ub_matrix((peak1, peak2))
    sample.ub_mat = ub_cal
    r_mat = tas.goniometer.r_mat(angles1)
    r_mat_cal = r_matrix_with_minimal_tilt(
        hkl=(0, 0, 2), ei=13.505137, ef=13.505137, two_theta=-51.530388, ub_mat=ub_matrix
    )

    angels1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2))
    for name, value in angels1_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-2)


def test_calculate_motor_angles():
    tas = TAS(fixed_ef=14.7)
    hb3_json = "./src/tavi/instrument/instrument_params/hb3_mnte.json"
    tas.load_instrument_params_from_json(hb3_json)
    lattice_params = (4.128474, 4.128474, 6.651507, 90, 90, 120)
    sample = Sample(lattice_params)
    sample.ub_mat = np.array(
        [
            [-0.000328, -0.004396, -0.150319],
            [-0.279308, -0.152332, 0.000314],
            [-0.014647, 0.234528, -0.002614],
        ]
    )
    # sample.plane_normal = [-0.017470, -0.052341, 0.998476]
    # sample.in_plane_ref = [-0.999847, 0.002086, -0.017385]
    tas.mount_sample(sample)

    angels1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2))
    angles1 = MotorAngles(two_theta=41.545383, omega=20.840750, sgl=1.001000, sgu=-3.000750)
    for name, value in angels1_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-2)

    # angels2_cal = tas.calculate_motor_angles(peak=(1.3, 0, 1.3), ei=ei + 35, ef=ef)
    # angles2 = MotorAngles(two_theta=21.02, omega=6.71, sgl=1.001000, sgu=-3.000750)
    # for name, value in angels2_cal._asdict().items():
    #     if value is not None:
    #         assert np.allclose(value, getattr(angles2, name), atol=1e-2)
