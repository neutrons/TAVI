import numpy as np
from tavi.utilities import *


source = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    #  Diameter D/4 if spherical
    "width": 15 / np.sqrt(12),
    "height": 15 / np.sqrt(12),
}


# guide before monochromator
guide = {
    "in_use": False,
    "div_h": 0.0,
    "div_v": 0.0,
}

monochromator = {
    "type": "PG002",
    "d_spacing": mono_ana_xtal["PG002"],
    "mosaic": 30 * min2rad,  # horizontal mosaic
    "mosaic_v": 30 * min2rad,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 7.62 / np.sqrt(12),
    "height": 10.16 / np.sqrt(12),
    "depth": 0.25 / np.sqrt(12),
    # horizontal focusing
    "curvh": 0.0,
    # vertical focusing
    "curvv": 0.0,
}

monitor = {
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 5 / np.sqrt(12),
    "height": 12 / np.sqrt(12),
}


analyzer = {
    "type": "Pg002",
    "d_spacing": mono_ana_xtal["Pg002"],
    "mosaic": 40 * min2rad,  # horizontal mosaic
    "mosaic_v": 40 * min2rad,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 7.62 / np.sqrt(12),
    "height": 7 / np.sqrt(12),
    "depth": 0.2 / np.sqrt(12),
    # horizontal focusing
    "curvh": 0.0,
    # vertical focusing
    "curvv": 0.0,
}
detector = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 4 / np.sqrt(12),
    "height": 12 / np.sqrt(12),
}

distances = {
    "src_mono": 650.0 * cm2A,
    "mono_sample": 190.0 * cm2A,
    "sample_ana": 160.0 * cm2A,
    "ana_det": 60.0 * cm2A,
    "mono_monitor": 86.0 * cm2A,
}

collimators = {  # in units of mins of arc
    "h_pre_mono": 48 * min2rad,
    "h_pre_sample": 40 * min2rad,
    "h_post_sample": 40 * min2rad,
    "h_post_ana": 120 * min2rad,
    "v_pre_mono": 150 * min2rad,
    "v_pre_sample": 270 * min2rad,
    "v_post_sample": 300 * min2rad,
    "v_post_ana": 600 * min2rad,
}


config_params = {
    "source": source,
    "guide": guide,
    "monochromator": monochromator,
    "monitor": monitor,
    "analyzer": analyzer,
    "detector": detector,
    "distances": distances,
    "collimators": collimators,
}
