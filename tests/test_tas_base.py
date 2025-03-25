import numpy as np

from tavi.instrument.tas import TAS
from tavi.sample import Sample


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
