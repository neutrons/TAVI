import matplotlib.pyplot as plt
import numpy as np
import pytest

from tavi.instrument.resolution.cooper_nathans import CN
from tavi.sample.xtal import Xtal

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, projection=None, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (3 * 2 * np.pi / tas.sample.c, 0, 0))
    assert rez.projection is None
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")


def test_hkl(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 0, 3))
    assert rez.projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert np.allclose(rez.angles, (60, 90, 90))
    assert rez.axes_labels == ("H (r.l.u.)", "K (r.l.u.)", "L (r.l.u.)", "E (meV)")


def test_projection(tas_params):
    tas, ei, ef, hkl, projection, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, projection=projection, R0=R0)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 3, 0))
    assert rez.projection == ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("(1, 1, 0)", "(0, 0, 1)", "(1, -1, 0)", "E (meV)")


def test_plot(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, R0=R0)
    rez.plot()

    plt.show()


@pytest.fixture
def tas_params():
    # cooper_nathans_CTAX

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = CN(SPICE_CONVENTION=False)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    nitio3 = Xtal.from_json(sample_json_path)
    tas.mount_sample(nitio3)

    ei = 4.8
    ef = 4.8
    hkl = (0, 0, 3)

    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = False

    tas_params = (tas, ei, ef, hkl, projection, R0)

    return tas_params
