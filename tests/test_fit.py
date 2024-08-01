import matplotlib.pyplot as plt
from tavi.tavi_data.tavi_data import TAVI_Data


def test_fit_scan(tavi):

    datasets = list(tavi.data.keys())[0]
    s = tavi.data[datasets]["scan0042"]
    (x, y, xerr, yerr, xlabel, ylabel, title, label) = s.generate_curve(
        norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25
    )

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    test_fit_scan(tavi)
