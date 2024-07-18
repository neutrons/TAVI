import numpy as np
import matplotlib.pylab as plt
from tavi.instrument.instrument_params.takin_test import instrument_params
from test_data_folder.test_samples.sample_test import test_xtal
from tavi.instrument.resolution.cooper_nathans import CN


def test_copper_nathans_localQ(tas):

    rez = tas.cooper_nathans(
        ei=13.5,
        ef=13.5,
        hkl=(0, 0, 2),
        projection=None,
        R0=False,
    )

    # describe and plot ellipses
    rez.calc_ellipses()
    rez.plot_ellipses()


def test_copper_nathans_hkl(tas):

    rez = tas.cooper_nathans(
        ei=13.5,
        ef=13.5,
        hkl=(0, 0, 2),
        R0=False,
    )

    # describe and plot ellipses
    rez.calc_ellipses()
    rez.plot_ellipses()


def test_copper_nathans_projection(tas):

    rez = tas.cooper_nathans(
        ei=13.5,
        ef=13.5,
        hkl=(0, 0, 2),
        projection=((1, 1, 0), (-1, 1, 0), (0, 0, 1)),
        R0=False,
    )

    # describe and plot ellipses
    rez.calc_ellipses()
    rez.plot_ellipses()


if __name__ == "__main__":

    tas = CN()
    tas.load_instrument(instrument_params)
    tas.load_sample(test_xtal)

    if tas.sample.ub_matrix is None:
        peak_list = [(0, 0, 2), (0, 2, 0)]
        angles_list = [
            (-51.530388, -45.220125, -0.000500, -2.501000),
            (-105.358735, 17.790125, -0.000500, -2.501000),
        ]
        tas.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)

    test_copper_nathans_localQ(tas)
    # test_copper_nathans_hkl(tas)
    # test_copper_nathans_projection(tas)

    plt.show()
