# from tavi.resolution.hb3 import config_params
from tavi.instrument.instrument_params.cg4c import cg4c_config_params
from tavi.instrument.instrument_params.takin_test import instrument_params
from test_data_folder.test_samples.nitio3 import nitio3
from test_data_folder.test_samples.sample_test import test_xtal
from tavi.utilities import *
from tavi.instrument.tas import TAS


def instrument_sample_setup(instrument_config, sample_config):

    tax = TAS()
    tax.load_instrument(instrument_config)
    tax.load_sample(sample_config)

    print(f"analyser type is {tax.analyzer.type}")
    print(f"a lattice parameter is {tax.sample.a}")
    print(f"gamma angle is {tax.sample.gamma}")

    return tax


def calc_ub_from_2_peaks_ctax(ctax):

    peak_list = [
        (0, 0, 3),
        (0.5, 0.5, 0),
    ]
    angles_list = [
        (53.240000, 32.865000, 2.310799, 2.021008),
        (48.489200, -59.507500, 2.310799, 2.021008),
    ]

    ctax.find_ub(peaks=peak_list, angles=angles_list, ei=4.799999, ef=4.799998)

    print(f"ub matrix = {ctax.sample.ub_matrix}")
    print(f"plane_normal = {ctax.sample.plane_normal}")
    print(f"in_plane_ref = {ctax.sample.in_plane_ref}")
    print(f"vec u = {ctax.sample.u}")
    print(f"vec v = {ctax.sample.v}")

    angles = ctax.find_angles(peak=(0, 0, 3), ei=4.799999, ef=4.799998)
    print(angles)
    angles = ctax.find_angles(peak=(0.5, 0.5, 0), ei=4.799999, ef=4.799998)
    print(angles)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)


def calc_ub_from_2_peaks_hb3(ctax):

    peak_list = [
        (0, 0, 2),
        (0, 2, 0),
    ]
    angles_list = [
        (-51.530388, -45.220125, -0.000500, -2.501000),
        (-105.358735, 17.790125, -0.000500, -2.501000),
    ]

    ctax.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)

    print(f"ub matrix = {np.round(ctax.sample.ub_matrix,6)}")
    print(f"plane_normal = {np.round(ctax.sample.plane_normal,6)}")
    print(f"in_plane_ref = {np.round(ctax.sample.in_plane_ref,6)}")
    print(f"vec u = {ctax.sample.u}")
    print(f"vec v = {ctax.sample.v}")


if __name__ == "__main__":

    # takin = instrument_sample_setup(instrument_params, test_xtal)
    # calc_ub_from_2_peaks_hb3(takin)

    ctax = instrument_sample_setup(cg4c_config_params, nitio3)
    calc_ub_from_2_peaks_ctax(ctax)
    pass
