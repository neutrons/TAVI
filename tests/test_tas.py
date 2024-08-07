from tavi.instrument.tas import TAS
from tavi.sample.xtal import Xtal
from tavi.utilities import *


def instrument_sample_setup(instrument_config, sample_config):
    tax = TAS()
    tax.load_instrument_from_dicts(instrument_config)
    tax.load_sample(sample_config)

    print(f"analyser type is {tax.analyzer.type}")
    print(f"a lattice parameter is {tax.sample.a}")
    print(f"gamma angle is {tax.sample.gamma}")

    return tax


def test_load_from_json(instrument_config_json_path, sample_json_path):
    # convert python dictionary to json file
    # with open("./src/tavi/instrument/instrument_params/takin_test.json", "w") as file:
    #     json.dump(instrument_config, file)

    # with open("./tests/test_data_folder/test_samples/nitio3.json", "w") as file:
    #     json.dump(nitio3, file)

    tax = TAS()
    tax.load_instrument_from_json(instrument_config_json_path)
    tax.load_sample_from_json(sample_json_path)

    print(f"analyser type is {tax.analyzer.type}")
    print(f"analyser type is {tax.analyzer.d_spacing}")
    print(f"a lattice parameter is {tax.sample.a}")
    print(f"gamma angle is {tax.sample.gamma}")

    return tax


def calc_ub_from_2_peaks_ctax(ctax):
    ub_matrix = np.array(
        [
            [-0.016934, -0.026164, -0.071871],
            [-0.108217, 0.12038, -0.003176],
            [0.20102, 0.192954, -0.007764],
        ]
    )
    plane_normal = [-0.04032, 0.998565, -0.035237]
    in_plane_ref = [-0.993257, -0.043892, -0.107299]

    peak_list = [
        (0, 0, 3),
        (0.5, 0.5, 0),
    ]
    # (s2, s1, sgl, sgu)
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


def calc_ub_from_2_peaks_hb3(hb3):
    lattice_params = (3.574924, 3.574924, 5.663212, 90, 90, 120)
    ub_matrix = np.array(
        [
            [0.053821, 0.107638, 0.166485],
            [0.272815, -0.013290, 0.002566],
            [0.164330, 0.304247, -0.058788],
        ]
    )
    b_matrix = np.array(
        [
            [0.3230, 0.1615, 0.0000],
            [0.0000, 0.2797, -0.0000],
            [0.0000, 0.0000, 0.1766],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    plane_normal = [0.000009, 0.999047, 0.043637]
    in_plane_ref = [0.942840, 0.014534, -0.332928]

    peak_list = [
        (0, 0, 2),
        (0, 2, 0),
    ]
    # (s2, s1, sgl, sgu)
    angles_list = [
        (-51.530388, -45.220125, -0.000500, -2.501000),
        (-105.358735, 17.790125, -0.000500, -2.501000),
    ]

    # hb3.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)
    hb3.find_ub(
        peaks=(peak_list[1], peak_list[0]),
        angles=(angles_list[1], angles_list[0]),
        ei=13.500172,
        ef=13.505137,
    )

    print(f"ub matrix = {np.round(hb3.sample.ub_matrix,6)}")
    print(f"plane_normal = {np.round(hb3.sample.plane_normal,6)}")
    print(f"in_plane_ref = {np.round(hb3.sample.in_plane_ref,6)}")
    print(f"vec u = {hb3.sample.u}")
    print(f"vec v = {hb3.sample.v}")

    angles = hb3.find_angles(peak=(0, 0, 2), ei=13.500172, ef=13.505137)
    print(angles)
    angles = hb3.find_angles(peak=(0, 2, 0), ei=13.500172, ef=13.505137)
    print(angles)


def calc_ub_from_2_peaks_hb1():
    hb1 = TAS()
    hb1.load_instrument_from_json("./src/tavi/instrument/instrument_params/hb3.json")

    lattice_params = (3.939520, 3.939520, 3.941957, 90.000000, 90.000000, 90.000000)
    xtal = Xtal(lattice_params=lattice_params)
    hb1.load_sample(xtal)
    ub_matrix = np.array(
        [
            [0.181283, 0.177659, -0.002718],
            [0.177656, -0.181301, -0.001380],
            [-0.002909, -0.000917, -0.253663],
        ]
    )
    b_matrix = np.array(
        [
            [0.253838, -0.000000, -0.000000],
            [0.000000, 0.253838, -0.000000],
            [0.000000, 0.000000, 0.253681],
        ]
    )

    plane_normal = [0.010097, 0.999933, -0.005550]
    in_plane_ref = [0.999892, -0.010155, -0.010658]

    peak_list = [
        (1, 1, 0),
        (0, 0, 1),
    ]
    # (s2, s1, sgl, sgu)
    angles_list = [
        (-52.449830, -26.825625, -0.578500, 0.318000),
        (-36.368334, -108.827250, -0.578500, 0.318000),
    ]

    # hb1.find_ub(peaks=peak_list, angles=angles_list, ei=13.499993, ef=13.506112)
    hb1.find_ub(
        peaks=(peak_list[1], peak_list[0]),
        angles=(angles_list[1], angles_list[0]),
        ei=13.499993,
        ef=13.506112,
    )
    print(f"ub matrix = {np.round(hb1.sample.ub_matrix,6)}")
    print(f"plane_normal = {np.round(hb1.sample.plane_normal,6)}")
    print(f"in_plane_ref = {np.round(hb1.sample.in_plane_ref,6)}")

    angles = hb1.find_angles(peak=(1, 1, 0), ei=13.499993, ef=13.506112)
    print(angles)
    angles = hb1.find_angles(peak=(0, 0, 1), ei=13.499993, ef=13.506112)
    print(angles)
    angles = hb1.find_angles(peak=(1, 1, 1), ei=13.499993, ef=13.506112)
    print(angles)


if __name__ == "__main__":
    takin = "./src/tavi/instrument/instrument_params/takin_test.json"
    sample_json_path = "./test_data/test_samples/nitio3.json"
    takin = test_load_from_json(takin, sample_json_path)
    calc_ub_from_2_peaks_hb3(takin)

    # ----------------
    cg4c = "./src/tavi/instrument/instrument_params/cg4c.json"
    nitio3 = "./test_data/test_samples/nitio3.json"
    ctax = test_load_from_json(cg4c, nitio3)
    calc_ub_from_2_peaks_ctax(ctax)

    # -----------------------------
    calc_ub_from_2_peaks_hb1()
