import numpy as np
import matplotlib.pylab as plt
from tavi.instrument.instrument_params.takin_test import instrument_params
from test_data_folder.test_samples.sample_test import test_xtal
from tavi.instrument.resolution.cooper_nathans import CN

from tavi.instrument.instrument_params.cg4c import cg4c_config_params
from test_data_folder.test_samples.nitio3 import nitio3

np.set_printoptions(floatmode="fixed", precision=4)


def test_copper_nathans_localQ(tas_params):

    tas, ei, ef, hkl, _, R0 = tas_params

    rez = tas.cooper_nathans(
        ei=ei,
        ef=ef,
        hkl=hkl,
        projection=None,
        R0=R0,
    )

    # describe and plot ellipses
    # rez.calc_ellipses()
    rez.plot()


def test_copper_nathans_hkl(tas_params):

    tas, ei, ef, hkl, _, R0 = tas_params

    rez = tas.cooper_nathans(
        ei=ei,
        ef=ef,
        hkl=hkl,
        R0=R0,
    )

    # describe and plot ellipses
    # rez.calc_ellipses()
    rez.plot()


def test_copper_nathans_projection(tas_params):

    tas, ei, ef, hkl, projection, R0 = tas_params

    rez = tas.cooper_nathans(
        ei=ei,
        ef=ef,
        hkl=hkl,
        projection=projection,
        R0=R0,
    )

    # describe and plot ellipses
    # rez.calc_ellipses()
    rez.plot()


def test_cooper_nathans_compare_3():

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

    ei = 13.5
    ef = 13.5
    hkl = (0, 0, 2)
    # projection = ((1, 0, 0), (-1, 2, 0), (0, 0, 1))
    projection = ((0, 0, 1), (0, -1, 0), (2, -1, 0))
    R0 = False

    tas_params = (tas, ei, ef, hkl, projection, R0)

    test_copper_nathans_localQ(tas_params)
    test_copper_nathans_hkl(tas_params)
    test_copper_nathans_projection(tas_params)

    plt.show()


def test_cooper_nathans_CTAX():

    tas = CN()
    tas.load_instrument(cg4c_config_params)
    tas.load_sample(nitio3)

    ei = 4.8
    ef = 4.8
    hkl = (0, 0, 3)

    projection = ((1, 1, 0), (0, 0, 1), (-1, 1, 0))
    R0 = False

    tas_params = (tas, ei, ef, hkl, projection, R0)

    test_copper_nathans_localQ(tas_params)
    test_copper_nathans_hkl(tas_params)
    test_copper_nathans_projection(tas_params)

    plt.show()


if __name__ == "__main__":

    test_cooper_nathans_compare_3()

    # test_cooper_nathans_CTAX()
