import numpy as np

from tavi.instrument.components.goni import Goniometer
from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import UBConf, plane_normal_from_two_peaks, uv_to_ub_matrix
from tavi.utilities import MotorAngles, Peak

if __name__ == "__main__":
    # cubic sample
    cube = Sample(lattice_params=(10, 10, 10, 90, 90, 90))
    u = (-1 / 2, np.sqrt(3) / 2, 0)
    v = (-np.sqrt(3) / 2, -1 / 2, 0)
    ub_mat = uv_to_ub_matrix(u, v, cube.lattice_params)
    u_mat = ub_mat.dot(np.linalg.inv(cube.b_mat))
    # expected values
    ub_mat_correct = (
        (-np.sqrt(3) / 20, -1 / 20, 0),
        (0, 0, 1 / 10),
        (-1 / 20, np.sqrt(3) / 20, 0),
    )
    b_mat_correct = (
        (1 / 10, 0, 0),
        (0, 1 / 10, 0),
        (0, 0, 1 / 10),
    )
    u_mat_correct = (
        (-np.sqrt(3) / 2, -1 / 2, 0),
        (0, 0, 1),
        (-1 / 2, np.sqrt(3) / 2, 0),
    )
    assert np.allclose(ub_mat, ub_mat_correct)
    assert np.allclose(cube.b_mat, b_mat_correct)
    assert np.allclose(u_mat, u_mat_correct)

    # create cg1b alignment station, using Mantid convention
    en = 14.0
    cg1b = TAS(fixed_ei=en, fixed_ef=en, convention=False)
    # conventioonal goniometer with s1 about +Y, sgl about+Z,  and sgu about -X
    # sense of s2 is negative
    cg1b.goniometer = Goniometer({"type": "Y,Z,-X", "sense": "-"})
    cg1b.mount_sample(cube)

    two_theta_100 = cg1b.get_two_theta(hkl=(-1, 0, 0))
    two_theta_100_correct = np.degrees(-2 * np.arcsin(9.045 / 2 / np.sqrt(en) / 10))
    assert np.allclose(two_theta_100, two_theta_100_correct, atol=1e-3)

    # check plane_normal and in_plane_ref defined by (1,0,0) and (0,1,0) peaks
    plane_normal, in_plnae_ref = plane_normal_from_two_peaks(u_mat, cube.b_mat, (1, 0, 0), (0, 1, 0))
    assert np.allclose(plane_normal, (0, 1, 0))
    assert np.allclose(in_plnae_ref, (-np.sqrt(3) / 2, 0, -1 / 2))

    cube.ub_conf = UBConf(ub_mat=ub_mat, plane_normal=plane_normal, in_plane_ref=in_plnae_ref)
    # check gonimoeter angles based on current UB matrix
    angles_m100 = cg1b.calculate_motor_angles(hkl=(-1, 0, 0))
    assert angles_m100 == MotorAngles(two_theta=-13.8845, omega=23.058, sgl=0, sgu=0)

    angles_010 = cg1b.calculate_motor_angles(hkl=(0, 1, 0))
    assert angles_010 == MotorAngles(two_theta=-13.8845, omega=113.058, sgl=0, sgu=0)

    ubconf_cal = cg1b.calculate_ub_matrix(
        peaks=(
            Peak((-1, 0, 0), angles_m100),
            Peak((0, 1, 0), angles_010),
        )
    )
    print(ubconf_cal)
    pass
