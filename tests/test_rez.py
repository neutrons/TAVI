import matplotlib.pylab as plt
from tavi.instrument.instrument_params.takin_test import instrument_params
from test_data_folder.test_samples.sample_test import test_xtal
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot1DManager, Plot2DManager
from tavi.utilities import *


def test_l_vs_en(tas):

    projection = ((0, 0, 1), (0, -1, 0), (2, -1, 0))

    q1 = [0, 2, 0.5]  # start, stop, step
    q2 = 0
    q3 = 0
    en = [0, 13, 2]  # start, stop, step

    ef = 13.5
    R0 = False

    plt2d = Plot2DManager()
    plt2d.rez_plot_gen(tas, projection, q1, q2, q3, en, ef, R0)


def test_h_vs_k(tas):

    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    q1 = [0, 2.1, 0.5]  # start, stop, step
    q2 = [0, 2.1, 0.5]  # start, stop, step
    q3 = 0
    en = 5
    ef = 13.5
    R0 = False

    plt2d = Plot2DManager()
    plt2d.rez_plot_gen(tas, projection, q1, q2, q3, en, ef, R0)


def test_1D(tas):

    projection = ((0, 0, 1), (0, -1, 0), (2, -1, 0))
    hkl = (0, 0, 1)
    en = 6

    ef = 13.5
    R0 = False

    axis = 1

    plt1d = Plot1DManager()
    plt1d.rez_plot_1D(tas, projection, hkl, en, ef, R0, axis)


if __name__ == "__main__":

    tas = CN()
    tas.load_instrument(instrument_params)
    tas.load_sample(test_xtal)

    # test_1D(tas)
    test_l_vs_en(tas)
    # test_h_vs_k(tas)

    plt.show()
