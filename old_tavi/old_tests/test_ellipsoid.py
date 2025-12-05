import matplotlib.pyplot as plt
import numpy as np
import pytest

from tavi.instrument.tas import TAS
from tavi.sample import Sample

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle, axes=None)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, 3 * 2 * np.pi / ctax.sample.c)
    assert rez.axes is None
    assert np.allclose(rez.angles, (90, 90, 90))


def test_hkl(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle=hkle)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q_vec, (0, 0, 3))
    assert rez.axes == ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")
    assert np.allclose(rez.angles, (60, 90, 90))


def test_projection(tas_params):
    ctax, hkle, projection = tas_params
    rez = ctax.cooper_nathans(hkle=hkle, axes=projection)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q_vec, (0, 3, 0))
    assert rez.axes == ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
    assert np.allclose(rez.angles, (90, 90, 90))


def test_plotting(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle=hkle)
    rez.plot()
    plt.show()


@pytest.fixture
def tas_params():
    # cooper_nathans_CTAX

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax = TAS(fixed_ef=4.8, convention="Spice")
    ctax.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    nitio3 = Sample.from_json(sample_json_path)
    ctax.mount_sample(nitio3)

    hkle = (0, 0, 3, 0)
    axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")

    tas_params = (ctax, hkle, axes)

    return tas_params
