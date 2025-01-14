import matplotlib.pyplot as plt
import numpy as np
import pytest
from mpl_toolkits.axisartist import Axes

from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot2D
from tavi.sample.xtal import Xtal

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(hkl_list=hkl, ei=ei, ef=ef, projection=None, R0=R0)
    ellipse = rez.get_ellipse(axes=(0, 3), PROJECTION=False)

    assert np.allclose(ellipse.angle, 90)
    assert ellipse.xlabel == "Q_para (1/A)"
    assert ellipse.ylabel == "E (meV)"


def test_hkl(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(hkl_list=hkl, ei=ei, ef=ef, R0=R0)

    e01_co = rez.get_ellipse(axes=(0, 1), PROJECTION=False)

    assert np.allclose(e01_co.angle, 60)
    assert e01_co.xlabel == "H (r.l.u.)"
    assert e01_co.ylabel == "K (r.l.u.)"


def test_plotting(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params
    rez = tas.cooper_nathans(hkl_list=hkl, ei=ei, ef=ef, R0=R0)

    e01_co = rez.get_ellipse(axes=(0, 1), PROJECTION=False)
    e01_inco = rez.get_ellipse(axes=(0, 1), PROJECTION=True)

    e03_co = rez.get_ellipse(axes=(0, 3), PROJECTION=False)
    e03_inco = rez.get_ellipse(axes=(0, 3), PROJECTION=True)

    p1 = Plot2D()
    p1.add_reso(e01_co, c="k", linestyle="solid")
    p1.add_reso(e01_inco, c="k", linestyle="dashed")

    p2 = Plot2D()
    p2.add_reso(e03_co, c="k", linestyle="solid", label="Coherent")
    p2.add_reso(e03_inco, c="k", linestyle="dashed", label="Incoherent")

    fig = plt.figure()
    ax1 = fig.add_subplot(
        121,
        axes_class=Axes,
        grid_helper=p1.grid_helper(e01_co.angle),
    )
    p1.plot(ax1)

    ax2 = fig.add_subplot(122, axes_class=Axes)
    p2.plot(ax2)

    fig.tight_layout(pad=2)
    plt.show()


@pytest.fixture
def tas_params():
    # cooper_nathans_CTAX

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = CooperNathans(SPICE_CONVENTION=False)
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
