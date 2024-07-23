import numpy as np

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# import scipy.constants as co
# ksq2E = (co.Planck / co.elementary_charge / 2.0 / np.pi) ** 2.0 * co.elementary_charge / 2.0 / co.neutron_mass * 1e23
ksq2eng = 2.072124855  # calculated with scipy.constants using the formula above

sig2fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
cm2angstrom = 1e8
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
def eng2k(en):
    """convert energy in meV to wave vector k in inverse Angstrom"""
    return


def get_angle(a, b, c):
    """In a triangle with sides a,b and c, get angle between a and b in radian
    Note:
        return value in [0,pi]"""
    acos = (a**2 + b**2 - c**2) / (2 * a * b)
    if acos > 1 or acos < -1:
        angle = None
    else:
        angle = np.arccos(acos)
    return angle


def get_angle_vec(v1, v2):
    """Get the angle in degress between two vectors v1 and v2"""
    return np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)) / np.pi * 180


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


def rot_x(nu):
    """rotation matrix about y-axis by angle nu

    Args:
        nu (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = nu / rad2deg
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [1, 0, 0],
            [0, c, -s],
            [0, s, c],
        ]
    )
    return mat


def rot_y(omega):
    """rotation matrix about y-axis by angle omega

    Args:
        omega (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = omega / rad2deg
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c],
        ]
    )
    return mat


def rot_z(mu):
    """rotation matrix about z-axis by angle mu

    Args:
        mu (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = mu / rad2deg
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1],
        ]
    )
    return mat
