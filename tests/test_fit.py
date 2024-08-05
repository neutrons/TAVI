import matplotlib.pyplot as plt
from tavi.tavi_data.tavi_data import TAVI_Data
from tavi.plotter import Plot1DManager


def test_fit_scan(tavi):

    datasets = list(tavi.data.keys())[0]
    s1 = tavi.data[datasets]["scan0042"]
    # x, y, xerr, yerr, xlabel, ylabel, title, label
    curve1 = s1.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    s2 = tavi.data[datasets]["scan0043"]
    curve2 = s2.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    p1 = Plot1DManager()
    p1.plot_curve(*curve1)
    p1.plot_curve(*curve2)

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    test_fit_scan(tavi)
