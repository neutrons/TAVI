import matplotlib.pyplot as plt
import numpy as np
import pytest

from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.tas import TAS
from tavi.sample import Sample

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle, projection=None)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (3 * 2 * np.pi / ctax.sample.c, 0, 0))
    assert rez.projection is None
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")


def test_hkl(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle=hkle)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 0, 3))
    assert rez.projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert np.allclose(rez.angles, (60, 90, 90))
    assert rez.axes_labels == ("(H, 0, 0) (r.l.u.)", "(0, K, 0) (r.l.u.)", "(0, 0, L) (r.l.u.)", "E (meV)")


def test_projection(tas_params):
    ctax, hkle, projection = tas_params
    rez = ctax.cooper_nathans(hkle=hkle, projection=projection)

    assert np.allclose(rez.hkl, (0, 0, 3))
    assert np.allclose(rez.q, (0, 3, 0))
    assert rez.projection == ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    assert np.allclose(rez.angles, (90, 90, 90))
    assert rez.axes_labels == ("(H, H, 0) (r.l.u.)", "(0, 0, K) (r.l.u.)", "(L, -L, 0) (r.l.u.)", "E (meV)")


def test_making_labels_from_projection():
    label = ResoEllipsoid.labels_from_projection()
    assert label == ("(H, 0, 0) (r.l.u.)", "(0, K, 0) (r.l.u.)", "(0, 0, L) (r.l.u.)", "E (meV)")

    label = ResoEllipsoid.labels_from_projection(projection=None)
    assert label == ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")

    label = ResoEllipsoid.labels_from_projection(projection=((1, 1, 0), (0, 0, 1), (1, -1, 0)))
    assert label == ("(H, H, 0) (r.l.u.)", "(0, 0, K) (r.l.u.)", "(L, -L, 0) (r.l.u.)", "E (meV)")

    label = ResoEllipsoid.labels_from_projection(projection=((1.0, 1.0, 0.0), (0, 0, 1), (1, -1, 0)))
    assert label == ("(H, H, 0) (r.l.u.)", "(0, 0, K) (r.l.u.)", "(L, -L, 0) (r.l.u.)", "E (meV)")


def test_plotting(tas_params):
    ctax, hkle, _ = tas_params
    rez = ctax.cooper_nathans(hkle=hkle)
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

    hkle = (0, 0, 3, 0)
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))

    tas_params = (ctax, hkle, projection)

    return tas_params
