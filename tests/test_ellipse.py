import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")  # temporarily disable interactive figure
import numpy as np
import pytest
from mpl_toolkits.axisartist import Axes

from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample

np.set_printoptions(floatmode="fixed", precision=4)


def test_local_q(tas_params):
    tas, hkl, _ = tas_params
    rez = tas.cooper_nathans(hkl=hkl, en=0, projection=None)
    ellipse = rez.get_ellipse(axes=(0, 3), PROJECTION=False)

    assert np.allclose(ellipse.angle, 90)
    assert ellipse.xlabel == "Q_para (1/A)"
    assert ellipse.ylabel == "E (meV)"


def test_hkl(tas_params):
    tas, hkl, _ = tas_params
    rez = tas.cooper_nathans(hkl=hkl, en=0)

    e01_co = rez.get_ellipse(axes=(0, 1), PROJECTION=False)

    assert np.allclose(e01_co.angle, 60)
    assert e01_co.xlabel == "H (r.l.u.)"
    assert e01_co.ylabel == "K (r.l.u.)"


def test_plotting(tas_params):
    tas, hkl, _ = tas_params
    rez = tas.cooper_nathans(hkl=hkl, en=0)

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
    tas = TAS(fixed_ef=4.8)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Sample.from_json(sample_json_path)
    tas.mount_sample(sample)

    hkl = (0, 0, 3)
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))

    tas_params = (tas, hkl, projection)

    return tas_params
