import numpy as np

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# import scipy.constants as co
# ksq2E = (co.Planck / co.elementary_charge / 2.0 / np.pi) ** 2.0 * co.elementary_charge / 2.0 / co.neutron_mass * 1e23
ksq2E = 2.072124855  # calculated with scipy.constants using the formula above

sig2fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
cm2A = 1e8
min2rad = 1.0 / 60.0 / 180.0 * np.pi
rad2deg = 180.0 / np.pi

# --------------------------------------------------------------------------
# d_spacing table from Shirane Appendix 3, in units of Angstrom, from
# --------------------------------------------------------------------------
mono_ana_xtal = {
    "PG002": 3.35416,
    "Pg002": 3.35416,
    "PG004": 1.67708,
    "Cu111": 2.08717,
    "Cu220": 1.27813,
    "Ge111": 3.26627,
    "Ge220": 2.00018,
    "Ge311": 1.70576,
    "Ge331": 1.29789,
    "Be002": 1.79160,
    "Be110": 1.14280,
    "Heusler": 3.435,  # Cu2MnAl(111)
}

# --------------------------------------------------------------------------
# helper functions
# --------------------------------------------------------------------------


def get_angle(v1, v2, v3):
    """return angle between v1 and v2 ni radian"""
    return np.arccos((v1**2 + v2**2 - v3**2) / (2 * v1 * v2))


def get_angle_bragg(q, d_spaceing):
    """return angle based on Bragg's law, in radian
    2d sin(theta) = lambda = 2 pi /q
    """
    return np.arcsin(np.pi / (d_spaceing * q))


def rotation_matrix_2d(phi):
    """rotate the coordination system by angle of phi about z-axis

    Args:
        phi (float): angle in radian

    Note:
        This is to rotate the coordination system, NOT the vector!


    """
    s = np.sin(phi)
    c = np.cos(phi)
    mat = np.array(
        [
            [c, s, 0],
            [-s, c, 0],
            [0, 0, 1],
        ]
    )
    return mat
