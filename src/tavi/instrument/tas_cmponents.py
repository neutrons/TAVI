import numpy as np
from tavi.utilities import *


class Source(object):
    """Neutron source"""

    def __init__(self, param_dict):

        self.name = None
        self.shape = None
        self.width = None
        self.height = None

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    setattr(self, key, val * cm2angstrom)
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
            if key in (
                "h_pre_mono",
                "h_pre_sample",
                "h_post_sample",
                "h_post_ana",
                "v_pre_mono",
                "v_pre_sample",
                "v_post_sample",
                "v_post_ana",
            ):
                setattr(self, key, val * min2rad)
            else:
                setattr(self, key, val)


# TODO
class Monitor(object):

    def __init__(self, param_dict):

        self.shape = "rectangular"
        self.width = 5
        self.height = 12

        for key, val in param_dict.items():
            setattr(self, key, val)

        if self.shape == "rectangular":
            self.width /= np.sqrt(12)
            self.height /= np.sqrt(12)
        elif self.shape == "spherical":
            pass


class Analyzer(object):

    def __init__(self, param_dict):

        self.type = "Pg002"
        self.d_spacing = mono_ana_xtal["Pg002"]
        self.mosaic = 45  # horizontal mosaic, min of arc
        self.mosaic_v = 45  # vertical mosaic, if anisotropic
        self.sense = -1
        # divide by np.sqrt(12) if rectangular
        # Diameter D/4 if spherical
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
                case "mosaic" | "mosaic_v" | "curvh" | "curvv":
                    setattr(self, key, val * min2rad)
                case "width" | "height" | "depth":
                    setattr(self, key, val * cm2angstrom)
                case _:
                    setattr(self, key, val)
            self.d_spacing = mono_ana_xtal[self.type]


class Detector(object):

    def __init__(self, param_dict):
        self.shape = "rectangular"
        self.width = 1.5  # in cm
        self.height = 5.0  # in cm

        for key, val in param_dict.items():
            setattr(self, key, val)

        # TODO
        if self.shape == "rectangular":
            self.width = self.width * cm2angstrom
            self.height = self.height * cm2angstrom
            # self.width /= np.sqrt(12)
            # self.height /= np.sqrt(12)
        elif self.shape == "spherical":
            pass  #
        else:
            print("unknown detector shape.")


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


class Monochromator(object):

    def __init__(self, param_dict):

        self.type = "PG002"
        self.d_spacing = mono_ana_xtal["PG002"]
        self.mosaic = 45  # horizontal mosaic, min of arc
        self.mosaic_v = 45  # vertical mosaic, if anisotropic
        self.sense = -1
        # divide by np.sqrt(12) if rectangular
        # Diameter D/4 if spherical
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
                case "mosaic" | "mosaic_v" | "curvh" | "curvv":
                    setattr(self, key, val * min2rad)
                case "width" | "height" | "depth":
                    setattr(self, key, val * cm2angstrom)
                case _:
                    setattr(self, key, val)
            self.d_spacing = mono_ana_xtal[self.type]

        for key, val in param_dict.items():
            setattr(self, key, val)


class Goniometer(object):
    """Goniometer table, type = TAS of 4C"""

    def __init__(self, param_dict):
        self.type = "TAS"
        self.sense = +1

        for key, val in param_dict.items():
            setattr(self, key, val)

    def r_mat(self, angles):
        "rotation matrix"
        # TODO
        if self.type == "TAS":  # Y-X-Z
            _, omega, _, _ = angles  # s1, s2, sgl. sgu
            r_mat = rot_y(omega)  # @ rot_x(-sgl) @ rot_z(-sgu)

        elif self.type == "4C":
            pass
        else:
            print("Unknow goniometer type. Needs to be TAS or 4C.")
        return r_mat

    def r_mat_inv(self, angles):
        """inverse of rotation matrix"""
        return np.linalg.inv(self.r_mat(angles))
