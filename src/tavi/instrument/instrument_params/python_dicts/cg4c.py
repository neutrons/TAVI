from tavi.utilities import *


source = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    #  Diameter D/4 if spherical
    "width": 7,  # / np.sqrt(12),  # in cm
    "height": 15,  # / np.sqrt(12),  # in cm
}


# guide before monochromator
# guide = {
#     "in_use": False,
#     "div_h": 0.0,
#     "div_v": 0.0,
# }

monochromator = {
    "type": "PG002",
    "d_spacing": mono_ana_xtal["PG002"],
    "mosaic": 30,  # horizontal mosaic
    "mosaic_v": 30,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 7.0,  # / np.sqrt(12),  # in cm
    "height": 15.0,  # / np.sqrt(12),  # in cm
    "depth": 0.2,  # / np.sqrt(12),  # in cm
    # horizontal focusing
    "curvh": 0.0,  # in cm^-1
    # vertical focusing
    "curvv": 60.4,  # at Ei = 4 meV, in cm^-1
}

monitor = {
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    # "width": 5 / np.sqrt(12),
    # "height": 12 / np.sqrt(12),
}

goniometer = {
    "sense": +1,
    "type": "Y-ZX",  # s1-sgl-sgu
    # CTAX's angle convention is actually YZ-X,
    # but the UB calculation in SPICE is done using the convetion Y-ZX
}

analyzer = {
    "type": "Pg002",
    "d_spacing": mono_ana_xtal["Pg002"],
    "mosaic": 90,  # horizontal mosaic
    "mosaic_v": 90,  # vertical mosaic, if anisotropic
    "sense": -1,
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 20.0 / np.sqrt(12),
    "height": 15.0 / np.sqrt(12),
    "depth": 0.2 / np.sqrt(12),
    # horizontal focusing
    "curvh": 0.0,
    # vertical focusing
    "curvv": 163.2,  # in cm^-1
}
detector = {
    "shape": "rectangular",  # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    "width": 5 / np.sqrt(12),
    "height": 10 / np.sqrt(12),
}

distances = {
    "src_mono": 530.0,  # in cm
    "mono_sample": 160.0,
    "sample_ana": 106.0,
    "ana_det": 50.0,
    # "mono_monitor": 86.0,
}

collimators = {  # in units of mins of arc
    "h_pre_mono": 40,
    "h_pre_sample": 100,
    "h_post_sample": 80,
    "h_post_ana": 120,
    "v_pre_mono": 600,
    "v_pre_sample": 600,
    "v_post_sample": 600,
    "v_post_ana": 600,
}


cg4c_config_params = {
    "source": source,
    # "guide": guide,
    "monochromator": monochromator,
    "goniometer": goniometer,
    "monitor": monitor,
    "analyzer": analyzer,
    "detector": detector,
    "distances": distances,
    "collimators": collimators,
}
