import matplotlib.pyplot as plt
import numpy as np
import pytest
from mpl_toolkits.axisartist import Subplot

from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot2D
from tavi.sample.xtal import Xtal

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, projection=None, R0=R0)
    ellipse = rez.get_ellipse(axes=(0, 3), PROJECTION=False)

    assert np.allclose(ellipse.angle, 90)
    assert ellipse.axes_labels == ("Q_para (1/A)", "E (meV)")


def test_hkl(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, R0=R0)
    ellipse = rez.get_ellipse(axes=(0, 1), PROJECTION=False)

    assert np.allclose(ellipse.angle, 60)
    assert ellipse.axes_labels == ("H (r.l.u.)", "K (r.l.u.)")

    p1 = Plot2D()
    p1.add_curve(ellipse.get_points(), fmt="-k")

    fig = plt.figure()
    ax = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax)
    p1.plot(ax)
    plt.show()


@pytest.fixture
def tas_params():
    # cooper_nathans_CTAX

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = CN(SPICE_CONVENTION=False)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Xtal.from_json(sample_json_path)
    tas.mount_sample(sample)

    ei = 4.8
    ef = 4.8
    hkl = (0, 0, 3)

    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = False

    tas_params = (tas, ei, ef, hkl, projection, R0)

    return tas_params
