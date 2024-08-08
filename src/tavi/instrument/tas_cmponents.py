import numpy as np

from tavi.utilities import *


class Source(object):
    """Neutron source"""

    def __init__(self, param_dict):
        self.name = None
        self.shape = "rectangular"
        self.width = None
        self.height = None

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    # divide by np.sqrt(12) if rectangular
                    #  Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "circular":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Source shape needs to be either rectangular or circular.")
                case _:
                    setattr(self, key, val)


class Guide(object):
    """Guide"""

    def __init__(self, param_dict):
        self.in_use = False
        self.div_h = None
        self.div_v = None

        for key, val in param_dict.items():
            match key:
                case "div_h" | "div_v":
                    setattr(self, key, val * min2rad)
                case _:
                    setattr(self, key, val)


class Monochromator(object):
    """Monochromator"""

    def __init__(self, param_dict):
        self.type = "PG002"
        self.d_spacing = mono_ana_xtal["PG002"]
        self.mosaic = 45  # horizontal mosaic, min of arc
        self.mosaic_v = 45  # vertical mosaic, if anisotropic
        self.sense = -1

        # divide by np.sqrt(12) if rectangular
        # Diameter D/4 if spherical
        self.shape = "rectangular"  # or spherical
        self.width = 12.0
        self.height = 8.0
        self.depth = 0.15
        # horizontal focusing
        self.curved_h = False
        self.curvh = 0.0
        self.optimally_curved_h = False
        # vertical focusing
        self.curved_v = False
        self.curvv = 0.0
        self.optimally_curved_v = False

        for key, val in param_dict.items():
            match key:
                # case "mosaic" | "mosaic_v" | "curvh" | "curvv":
                #     setattr(self, key, val * min2rad)
                case "width" | "height" | "depth":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Monochromator shape needs to be either rectangular or spherical.")
                case _:
                    setattr(self, key, val)
            self.d_spacing = mono_ana_xtal[self.type]


class Collimators(object):
    """collimitor divergences, in mins of arc"""

    def __init__(self, param_dict):
        self.h_pre_mono = 30  # mins of arc
        self.h_pre_sample = 30
        self.h_post_sample = 30
        self.h_post_ana = 30
        self.v_pre_mono = 30
        self.v_pre_sample = 30
        self.v_post_sample = 30
        self.v_post_ana = 30

        for key, val in param_dict.items():
            # if key in (
            #     "h_pre_mono",
            #     "h_pre_sample",
            #     "h_post_sample",
            #     "h_post_ana",
            #     "v_pre_mono",
            #     "v_pre_sample",
            #     "v_post_sample",
            #     "v_post_ana",
            # ):
            #     setattr(self, key, val * min2rad)
            # else:
            setattr(self, key, val)


# TODO
class Monitor(object):
    def __init__(self, param_dict):
        self.shape = "rectangular"  # or spherical
        self.width = 5 * cm2angstrom
        self.height = 12 * cm2angstrom

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Monitor shape needs to be either rectangular or spherical.")
                case _:
                    setattr(self, key, val)


class Analyzer(object):
    def __init__(self, param_dict):
        self.type = "Pg002"
        self.d_spacing = mono_ana_xtal["Pg002"]
        self.mosaic = 45  # horizontal mosaic, min of arc
        self.mosaic_v = 45  # vertical mosaic, if anisotropic
        self.sense = -1
        # divide by np.sqrt(12) if rectangular
        # Diameter D/4 if spherical
        self.shape = "rectangular"
        self.width = 12.0
        self.height = 8.0
        self.depth = 0.3
        # horizontal focusing
        self.curved_h = False
        self.curvh = 0.0
        self.optimally_curved_h = False
        # vertical focusing
        self.curved_v = False
        self.curvv = 0.0
        self.optimally_curved_v = False

        for key, val in param_dict.items():
            match key:
                # case "mosaic" | "mosaic_v" | "curvh" | "curvv":
                #     setattr(self, key, val * min2rad)
                case "width" | "height" | "depth":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Analyzer shape needs to be either rectangular or spherical.")

                    setattr(self, key, val * cm2angstrom)
                case _:
                    setattr(self, key, val)
            self.d_spacing = mono_ana_xtal[self.type]


# TODO
class Detector(object):
    def __init__(self, param_dict):
        self.shape = "rectangular"
        self.width = 1.5  # in cm
        self.height = 5.0  # in cm

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Detector shape needs to be either rectangular or spherical.")
                case _:
                    setattr(self, key, val)


class Arms(object):
    """lengths of arms, in units of cm"""

    def __init__(self, param_dict):
        # in units of cm
        self.src_mono = 0
        self.mono_sample = 0
        self.sample_ana = 0
        self.ana_det = 0
        self.mono_monitor = 0

        for key, val in param_dict.items():
            setattr(self, key, val * cm2angstrom)


class Goniometer(object):
    """Goniometer table, type = Y-ZX or YZ-X"""

    def __init__(self, param_dict):
        self.type = "Y-ZX"  # Y-mZ-X for Huber stage at HB1A and HB3
        self.sense = -1

        for key, val in param_dict.items():
            setattr(self, key, val)

    def r_mat(self, angles):
        "rotation matrix"

        omega, sgl, sgu = angles  # s2, s1, sgl, sgu
        match self.type:
            case "Y-ZX":  # HB3
                r_mat = rot_y(omega) @ rot_z(-1 * sgl) @ rot_x(sgu)
            case "YZ-X":  # CG4C
                r_mat = rot_y(omega) @ rot_z(sgl) @ rot_x(-1 * sgu)
            case _:
                r_mat = None
                print("Unknow goniometer type. Curruntly support Y-ZX and YZ-X")

        return r_mat

    def r_mat_inv(
        self,
        angles,
    ):
        """inverse of rotation matrix"""
        # return np.linalg.inv(self.r_mat(angles))
        return self.r_mat(angles).T

    def angles_from_r_mat(self, r_mat):
        """Calculate goniometer angles from the R matrix

        Note:
            range of np.arcsin is -pi/2 to pi/2
            range of np.atan2 is -pi to pi
        """

        match self.type:
            case "Y-ZX" | "YZ-X":  # Y-mZ-X (s1, sgl, sgu) for HB1A and HB3, Y-Z-mX (s1, sgl, sgu) for CG4C
                # sgl1 = np.arcsin(r_mat[1, 0]) * rad2deg
                # sgl2 = np.arccos(np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
                sgl = np.arctan2(r_mat[1, 0], np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg

                # sgu1 = np.arcsin(-r_mat[1, 2] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
                # sgu2 = np.arccos(r_mat[1, 1] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
                sgu = (
                    np.arctan2(
                        -r_mat[1, 2] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                        r_mat[1, 1] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                    )
                    * rad2deg
                )

                # omega1 = np.arcsin(-r_mat[2, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
                # omega2 = np.arccos(r_mat[0, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
                omega = (
                    np.arctan2(
                        -r_mat[2, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                        r_mat[0, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                    )
                    * rad2deg
                )
                match self.type:
                    case "Y-ZX":
                        angles = (omega, -1 * sgl, sgu)
                    case "YZ-X":
                        angles = (omega, sgl, -1 * sgu)

            case _:
                angles = None
                print("Unknow goniometer type.  Curruntly support Y-ZX and YZ-X.")

        return angles
