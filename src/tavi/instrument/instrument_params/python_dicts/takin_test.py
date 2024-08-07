from tavi.utilities import *

source = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    #  Diameter D/4 if spherical
    "width": 6.0,  # in cm
    "height": 12.0,  # in cm
}


# guide before monochromator
guide = {
    "in_use": False,
    "div_h": 15.0,  # min of arc
    "div_v": 15.0,  # min of arc
}

monochromator = {
    "type": "PG002",
    "shape": "rectangular",
    # "d_spacing": mono_ana_xtal["PG002"],
    "mosaic": 45,  # horizontal mosaic
    "mosaic_v": 45,  # vertical mosaic, if anisotropic
    "sense": +1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 12.0,
    "height": 8.0,
    "depth": 0.15,
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
    "shape": "rectangular",
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 5,  # / np.sqrt(12),
    "height": 12,  # / np.sqrt(12),
}

goniometer = {
    "sense": -1,
    "type": "Y-ZX",  # s1-sgl-sgu
}

analyzer = {
    "type": "Pg002",
    "d_spacing": mono_ana_xtal["Pg002"],
    "mosaic": 45,  # horizontal mosaic
    "mosaic_v": 45,  # vertical mosaic, if anisotropic
    "sense": +1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 12.0,
    "height": 8.0,
    "depth": 0.3,
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
    "width": 1.5,  # cm
    "height": 5.0,  # cm
}

distances = {  # in units of cm
    "src_mono": 10.0,
    "mono_sample": 200.0,
    "sample_ana": 115.0,
    "ana_det": 85.0,
    # "mono_monitor": 86.0 ,
}

collimators = {  # in units of mins of arc
    "h_pre_mono": 30,
    "h_pre_sample": 30,
    "h_post_sample": 30,
    "h_post_ana": 30,
    "v_pre_mono": 30,
    "v_pre_sample": 30,
    "v_post_sample": 30,
    "v_post_ana": 30,
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
