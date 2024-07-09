import numpy as np
import matplotlib.pyplot as plt
from tavi.tavi_data.spice_to_nexus import convert_spice_to_nexus
from tavi.tavi_data.tavi_data import TAVI_Data


def test_conversion(exp_numer):

    spice_folder = f"./tests/test_data_folder/exp{exp_numer}/"
    # h5_file_name = "./tests/test_data_folder/tavi_exp758.h5"
    nexus_file_name = f"./tests/test_data_folder/nexus_exp{exp_numer}.h5"
    convert_spice_to_nexus(spice_folder, nexus_file_name)


def test_load_nexus_to_new_tavi(tavi):

    tavi_file_name = "./tests/test_data_folder/tavi_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_file_name = "./tests/test_data_folder/nexus_exp424.h5"
    tavi.load_tavi_data_from_disk(nexus_file_name)
    return tavi


def test_plot_scan(tavi):

    print(len(tavi.data))
    s = tavi.data["scan0025"]
    x, y, xerr, yerr, xlabel, ylabel, title = s.curve_gen()

    fig, ax = plt.subplots()
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.6)

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    # test_conversion(424)
    test_load_nexus_to_new_tavi(tavi)

    test_plot_scan(tavi)
