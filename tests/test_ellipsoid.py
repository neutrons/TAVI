import matplotlib.pyplot as plt
import numpy as np
import pytest

from tavi.instrument.tas import TAS
from tavi.sample import Sample

np.set_printoptions(floatmode="fixed", precision=4)


def test_out_of_reach(tas_params):
    ctax, hkl, _, R0 = tas_params
    rez = ctax.cooper_nathans(hkl=(0, 0, 0), en=0, projection=None, R0=R0)
    assert np.allclose(rez.STATUS, False)

    rez = ctax.cooper_nathans(hkl=(10, 10, 10), en=0, projection=None, R0=R0)
    assert np.allclose(rez.STATUS, False)


def test_local_q(tas_params):
    ctax, hkl, _, R0 = tas_params
    rez = ctax.cooper_nathans(hkl, en=0, projection=None, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (3 * 2 * np.pi / ctax.sample.c, 0, 0))
    assert rez.projection is None
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")


def test_hkl(tas_params):
    ctax, hkl, _, R0 = tas_params
    rez = ctax.cooper_nathans(hkl=hkl, en=0, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 0, 3))
    assert rez.projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert np.allclose(rez.angles, (60, 90, 90))
    assert rez.axes_labels == ("H (r.l.u.)", "K (r.l.u.)", "L (r.l.u.)", "E (meV)")


def test_projection(tas_params):
    ctax, hkl, projection, R0 = tas_params
    rez = ctax.cooper_nathans(hkl=hkl, en=0, projection=projection, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 3, 0))
    assert rez.projection == ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("(1, 1, 0)", "(0, 0, 1)", "(1, -1, 0)", "E (meV)")


def test_plotting(tas_params):
    ctax, hkl, _, R0 = tas_params
    rez = ctax.cooper_nathans(hkl=hkl, en=0, R0=R0)
    rez.plot_ellipses()
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

    hkl = (0, 0, 3)
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = False

    tas_params = (ctax, hkl, projection, R0)

    return tas_params
