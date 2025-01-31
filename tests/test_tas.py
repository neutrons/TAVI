import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample.xtal import Xtal
from tavi.utilities import MotorAngles


def test_find_two_theta():
    ctax = TAS()
    ctax_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(ctax_json)
    nitio3 = Xtal.from_json("./test_data/test_samples/nitio3.json")
    ctax.mount_sample(nitio3)

    two_theta = ctax.get_two_theta(hkl=(0, 0, 3), ei=4.799999)
    assert np.allclose(two_theta, 53.240000, atol=1e-1)
    two_theta = ctax.get_two_theta(hkl=(0.5, 0.5, 0), ei=4.799999)
    assert np.allclose(two_theta, 48.489200, atol=1e-1)


def test_calculate_motor_angles():
    tas = TAS()
    hb3_json = "./src/tavi/instrument/instrument_params/hb3_mnte.json"
    tas.load_instrument_params_from_json(hb3_json)
    lattice_params = (4.128474, 4.128474, 6.651507, 90, 90, 120)
    sample = Xtal(lattice_params)
    sample.ub_mat = np.array(
        [
            [-0.000328, -0.004396, -0.150319],
            [-0.279308, -0.152332, 0.000314],
            [-0.014647, 0.234528, -0.002614],
        ]
    )
    sample.plane_normal = [-0.017470, -0.052341, 0.998476]
    sample.in_plane_ref = [-0.999847, 0.002086, -0.017385]
    tas.mount_sample(sample)

    ei = 14.7
    ef = 14.7

    angels1_cal = tas.calculate_motor_angles(peak=(0, 0, 2), ei=ei, ef=ef)
    angles1 = MotorAngles(two_theta=41.545383, omega=20.840750, sgl=1.001000, sgu=-3.000750)
    for name, value in angels1_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles1, name), atol=1e-2)

    angels2_cal = tas.calculate_motor_angles(peak=(1.3, 0, 1.3), ei=ei + 35, ef=ef)
    angles2 = MotorAngles(two_theta=21.02, omega=6.71, sgl=1.001000, sgu=-3.000750)
    for name, value in angels2_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles2, name), atol=1e-2)
