import numpy as np
import pytest

from tavi.instrument.tas import TAS
from tavi.sample import Sample


@pytest.fixture
def param_dict():
    return {
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


def test_load_from_dict(param_dict):
    tax = TAS()
    tax._load_instrument_parameters(param_dict)

    assert tax.analyzer.type == "Pg002"


def test_load_from_json():
    tax = TAS()
    cg4c = "./src/tavi/instrument/instrument_params/cg4c.json"
    nitio3 = Sample.from_json("./test_data/test_samples/nitio3.json")
    tax.load_instrument_params_from_json(cg4c)
    tax.mount_sample(nitio3)

    assert tax.analyzer.type == "Pg002"
    assert np.allclose(tax.analyzer.d_spacing, 3.35416)
    assert np.allclose(tax.sample.a, 5.034785)
    assert np.allclose(tax.sample.gamma, 120)


def test_set_motor_range():
    ctax = TAS()
    cg4c_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax.load_instrument_params_from_json(cg4c_json)
    ctax.goniometer.set_limit(motor_name="omega", motor_range=(-200, 200))
    assert ctax.goniometer.limits["omega"] == (-200, 200)
