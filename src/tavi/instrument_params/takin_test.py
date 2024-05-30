import numpy as np
from tavi.utilities import *


source = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    #  Diameter D/4 if spherical
    "width": 6.0 * cm2A,
    "height": 12.0 * cm2A,
}


# guide before monochromator
guide = {
    "in_use": False,
    "div_h": 15.0 * min2rad,
    "div_v": 15.0 * min2rad,
}

monochromator = {
    "type": "PG002",
    "d_spacing": mono_ana_xtal["PG002"],
    "mosaic": 45 * min2rad,  # horizontal mosaic
    "mosaic_v": 45 * min2rad,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 12.0 * cm2A,
    "height": 8.0 * cm2A,
    "depth": 0.15 * cm2A,
    # horizontal focusing
    "curved_h": False,
    "curvh": 0.0,
    "optimally_curved_h": False,
    # vertical focusing
    "curved_v": False,
    "curvv": 0.0,
    "optimally_curved_v": False,
}

monitor = {
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 5 / np.sqrt(12),
    "height": 12 / np.sqrt(12),
}

goniometer = {
    "sense": +1,
    "num_axis": "TAS",  # "TAS" or "4C"
}

analyzer = {
    "type": "Pg002",
    "d_spacing": mono_ana_xtal["Pg002"],
    "mosaic": 45 * min2rad,  # horizontal mosaic
    "mosaic_v": 45 * min2rad,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 12.0 * cm2A,
    "height": 8.0 * cm2A,
    "depth": 0.3 * cm2A,
    # horizontal focusing
    "curved_h": False,
    "curvh": 0.0,
    "optimally_curved_h": False,
    # vertical focusing
    "curved_v": False,
    "curvv": 0.0,
    "optimally_curved_v": False,
}
detector = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 1.5 * cm2A,
    "height": 5.0 * cm2A,
}

distances = {
    "src_mono": 10.0 * cm2A,
    "mono_sample": 200.0 * cm2A,
    "sample_ana": 115.0 * cm2A,
    "ana_det": 85.0 * cm2A,
    # "mono_monitor": 86.0 * cm2A,
}

collimators = {  # in units of mins of arc
    "h_pre_mono": 30 * min2rad,
    "h_pre_sample": 30 * min2rad,
    "h_post_sample": 30 * min2rad,
    "h_post_ana": 30 * min2rad,
    "v_pre_mono": 30 * min2rad,
    "v_pre_sample": 30 * min2rad,
    "v_post_sample": 30 * min2rad,
    "v_post_ana": 30 * min2rad,
}


instrument_params = {
    "source": source,
    "guide": guide,
    "monochromator": monochromator,
    "monitor": monitor,
    "goniometer": goniometer,
    "analyzer": analyzer,
    "detector": detector,
    "distances": distances,
    "collimators": collimators,
}
