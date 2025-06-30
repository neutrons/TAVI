import copy

import numpy as np
import pytest

from tavi.instrument.resolution.resolution_calculator import ResolutionCalculator
from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import UBConf


@pytest.fixture
def params():
    instrument_params = {
        "source": {"shape": "rectangular", "width": 7.0, "height": 15.0},
        "guide": {"in_use": False, "div_h": 0.0, "div_v": 0.0},
        "monochromator": {
            "type": "PG002",
            "mosaic_h": 30,
            "mosaic_v": 30,
            "sense": "-",
            "shape": "rectangular",
            "width": 7,
            "height": 15,
            "depth": 0.2,
            "curved_h": False,
            "curvh": 0.0,
            "optimally_curved_h": False,
            "curved_v": True,
            "curvv": 60.4,
            "optimally_curved_v": False,
        },
        "goniometer": {"sense": "+", "type": "Y,-Z,X"},
        "monitor": {"shape": "rectangular"},
        "analyzer": {
            "type": "Pg002",
            "d_spacing": 3.35416,
            "mosaic_h": 90,
            "mosaic_v": 90,
            "sense": "-",
            "shape": "rectangular",
            "width": 20.0,
            "height": 18.0,
            "depth": 0.2,
            "curved_h": False,
            "curvh": 0.0,
            "optimally_curved_h": False,
            "curved_v": True,
            "curvv": 163.2,
            "optimally_curved_v": False,
        },
        "detector": {"shape": "rectangular", "width": 5, "height": 10.0},
        "distances": {"src_mono": 530.0, "mono_sample": 160.0, "sample_ana": 106.0, "ana_det": 50.0},
        "collimators": {
            "h_pre_mono": 40,
            "h_pre_sample": 100,
            "h_post_sample": 80,
            "h_post_ana": 120,
            "v_pre_mono": 600,
            "v_pre_sample": 600,
            "v_post_sample": 600,
            "v_post_ana": 600,
        },
    }

    sample = Sample((5.034785, 5.034785, 13.812004, 90, 90, 120))
    sample.set_mosaic(30, 30)
    sample.ub_conf = UBConf(
        ub_mat=np.array(
            [
                [-0.016965, -0.026212, -0.071913],
                [-0.201388, -0.193307, 0.007769],
                [-0.108415, 0.120600, -0.003178],
            ]
        ),
        plane_normal=[-0.04032, 0.035237, 0.998565],
        in_plane_ref=[-0.993257, 0.107299, -0.043892],
    )

    return instrument_params, sample


@pytest.fixture
def ctax():
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax = TAS(fixed_ef=4.8)
    ctax.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Sample.from_json(sample_json_path)
    ctax.mount_sample(sample)

    return ctax


def test_validate_instrument_parameter(ctax, params):
    rc = ResolutionCalculator(instrument=ctax)
    rc.validate_instrument_parameters()

    instrument_params, sample = params
    tas = TAS()
    tas.mount_sample(sample)
    with pytest.raises(AttributeError) as e_info:
        rc = ResolutionCalculator(instrument=tas)
        rc.validate_instrument_parameters()
    assert "monochromator is missing in TAS(fixed_ei=None, fixed_ef=None, convention=Spice)." in str(e_info.value)

    p1 = copy.deepcopy(instrument_params)
    p1["monochromator"].pop("mosaic_h")
    tas._load_instrument_parameters(p1)
    with pytest.raises(AttributeError) as e_info:
        rc = ResolutionCalculator(instrument=tas)
        rc.validate_instrument_parameters()
    assert "mosaic_h is missing in Monochromator/Analyzer of type PG002." in str(e_info.value)

    p2 = copy.deepcopy(instrument_params)
    p2["monochromator"]["mosaic_h"] = -1
    tas._load_instrument_parameters(p2)
    with pytest.raises(ValueError) as e_info:
        rc = ResolutionCalculator(instrument=tas)
        rc.validate_instrument_parameters()
    assert "mosaic_h = -1 in Monochromator/Analyzer of type PG002 cannot be negative." in str(e_info.value)

    p3 = copy.deepcopy(instrument_params)
    p3["collimators"]["h_pre_mono"] = -1
    tas._load_instrument_parameters(p3)
    with pytest.raises(ValueError) as e_info:
        rc = ResolutionCalculator(instrument=tas)
        rc.validate_instrument_parameters()
    assert "horizontal_divergence = [-1, 100, 80, 120] in Collimator cannot be negative." in str(e_info.value)

    p4 = copy.deepcopy(instrument_params)
    p4["goniometer"].pop("sense")
    tas._load_instrument_parameters(p4)
    with pytest.raises(AttributeError) as e_info:
        rc = ResolutionCalculator(instrument=tas)
        rc.validate_instrument_parameters()
    assert "sense is missing in Goniometer, type=Y,-Z,X." in str(e_info.value)


def test_generate_hkle_list():
    hkle = ResolutionCalculator.generate_hkle(
        axes=((1, 1, 0), (0, 0, 1), (1, -1, 0), "en"),
        grid=((-0.5, 0.15, 0.05), 3, 0, (0, 4.1, 0.4)),
    )
    assert len(hkle) == 143
