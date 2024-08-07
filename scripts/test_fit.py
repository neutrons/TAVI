import matplotlib.pyplot as plt
from tavi.tavi_data.tavi_data import TAVI_Data
from tavi.tavi_data.fit import Fit
from tavi.plotter import Plot1DManager


def test_fit_scan(tavi):

    p1 = Plot1DManager()

    datasets = list(tavi.data.keys())[0]
    s1 = tavi.data[datasets]["scan0042"]
    curve1 = s1.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    x, y, xerr, yerr, xlabel, ylabel, title, label = curve1
    f1 = Fit(x=x, y=y, err=yerr, fit_range=(0.5, 3))
    f1.add_background()
    f1.add_signal()
    # f1.add_signal()
    f1.perform_fit()

    p1.plot_curve(*curve1)

    # s2 = tavi.data[datasets]["scan0043"]
    # curve2 = s2.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)
    # p1.plot_curve(*curve2)

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    test_fit_scan(tavi)
