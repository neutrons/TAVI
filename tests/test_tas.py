import numpy as np

from tavi.instrument.components.goni import Goniometer
from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import plane_normal_from_two_peaks, r_matrix_with_minimal_tilt, uv_to_ub_matrix
from tavi.utilities import MotorAngles, Peak, UBConf, mantid_to_spice, spice_to_mantid


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

    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

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
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

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
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

    angles_cal = tas.goniometer.angles_from_r_mat(r_mat, two_theta=angles.two_theta)
    for name, value in angles_cal._asdict().items():
        if value is not None:
            assert np.allclose(value, getattr(angles, name), atol=1e-2)


def test_calculate_ub():

    tas = TAS(fixed_ef=13.505137, spice_convention=False)
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

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
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

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
    tas.goniometer = Goniometer({"sense": "-", "type": "Y,-Z,X"})

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
    assert np.allclose(r_mat, spice_to_mantid(r_mat_cal.T).T, atol=1e-3)

    angles1_cal = tas.calculate_motor_angles(hkl=(0, 0, 2))
    assert angles1_cal == angles1

    angles2_cal = tas.calculate_motor_angles(hkl=(1, 0, 0))
    angles2 = MotorAngles(two_theta=38.526091, omega=-70.614125, sgl=1.001000, sgu=-3.000750)
    assert angles2_cal == angles2

    angles3_cal = tas.calculate_motor_angles(hkl=(1.3, 0, 1.3))
    assert np.allclose(angles3_cal.sgl, 1.001000, atol=1e-2)
    assert np.allclose(angles3_cal.sgu, -3.000750, atol=1e-2)


def test_cube():
    a, b, c = 10, 20, 30
    cube = Sample(lattice_params=(a, b, c, 90, 90, 90))
    u = (1, 0, 0)
    v = (0, 1, 0)
    ub_mat = uv_to_ub_matrix(u, v, cube.lattice_params)
    assert np.allclose(ub_mat, ((0, 1 / b, 0), (0, 0, 1 / c), (1 / a, 0, 0)))
    assert np.allclose(cube.b_mat, ((1 / a, 0, 0), (0, 1 / b, 0), (0, 0, 1 / c)))
    u_mat = ub_mat.dot(np.linalg.inv(cube.b_mat))
    assert np.allclose(u_mat, ((0, 1, 0), (0, 0, 1), (1, 0, 0)))

    hb1a = TAS(fixed_ei=14.45, fixed_ef=14.45)
    hb1a.mount_sample(cube)

    # TAS goniometer with sgl and sgu
    hb1a.goniometer = Goniometer({"type": "Y,-Z,X", "sense": "-"})
    two_theta_100 = hb1a.get_two_theta(hkl=(1, 0, 0))
    assert np.allclose(two_theta_100, np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / a)), atol=1e-3)
    two_theta_010 = hb1a.get_two_theta(hkl=(0, 1, 0))
    assert np.allclose(two_theta_010, np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / b)), atol=1e-3)
    two_theta_110 = hb1a.get_two_theta(hkl=(1, 1, 0))
    assert np.allclose(
        two_theta_110,
        np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(14.45) / (a * b) * np.sqrt(a**2 + b**2))),
        atol=1e-3,
    )

    plane_normal, in_plnae_ref = plane_normal_from_two_peaks(u_mat, cube.b_mat, (1, 0, 0), (0, 1, 0))
    cube.ub_conf = UBConf(
        ub_mat=mantid_to_spice(ub_mat),
        plane_normal=mantid_to_spice(plane_normal),
        in_plane_ref=mantid_to_spice(in_plnae_ref),
    )

    angles_100_tas = hb1a.calculate_motor_angles(hkl=(1, 0, 0))
    assert np.allclose(angles_100_tas.two_theta, two_theta_100)
    assert np.allclose(angles_100_tas.omega, hb1a.get_psi((1, 0, 0)))
    assert np.allclose(angles_100_tas.sgl, 0)
    assert np.allclose(angles_100_tas.sgu, 0)

    angles_010_tas = hb1a.calculate_motor_angles(hkl=(0, 1, 0))
    assert np.allclose(angles_010_tas.two_theta, two_theta_010)
    assert np.allclose(angles_010_tas.omega, -90 + hb1a.get_psi((0, 1, 0)))
    assert np.allclose(angles_010_tas.sgl, 0)
    assert np.allclose(angles_010_tas.sgu, 0)

    angles_110_tas = hb1a.calculate_motor_angles(hkl=(1, 1, 0))
    assert np.allclose(angles_110_tas.two_theta, two_theta_110)
    assert np.allclose(
        angles_110_tas.omega,
        np.degrees(-np.arctan(1 / 2)) + hb1a.get_psi((1, 1, 0)),
    )
    assert np.allclose(angles_110_tas.sgl, 0)
    assert np.allclose(angles_110_tas.sgu, 0)

    # 4C goniometer with chi and phi
    hb1a.goniometer = Goniometer({"type": "Y,Z,Y,bisect", "sense": "-"})
    angles_100_4c = hb1a.calculate_motor_angles(hkl=(1, 0, 0))
    assert np.allclose(angles_100_4c.two_theta, two_theta_100)
    assert np.allclose(angles_100_4c.omega, two_theta_100 / 2)
    assert np.allclose(angles_100_4c.chi, 0)
    assert np.allclose(angles_100_4c.phi, 90)

    angles_010_4c = hb1a.calculate_motor_angles(hkl=(0, 1, 0))
    assert np.allclose(angles_010_4c.two_theta, two_theta_010)
    assert np.allclose(angles_010_4c.omega, two_theta_010 / 2)
    assert np.allclose(angles_010_4c.chi, 0)
    assert np.allclose(angles_010_4c.phi, 0)

    angles_110_4c = hb1a.calculate_motor_angles(hkl=(1, 1, 0))
    assert np.allclose(angles_110_4c.two_theta, two_theta_110)
    assert np.allclose(angles_110_4c.omega, two_theta_110 / 2)
    assert np.allclose(angles_110_4c.chi, 0)
    assert np.allclose(
        angles_110_4c.phi,
        np.degrees(-np.arctan(1 / 2)) + hb1a.get_psi((1, 1, 0)) - two_theta_110 / 2,
    )
