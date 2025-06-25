import matplotlib.pyplot as plt
import numpy as np

from tavi.data.scan import Scan


def test_scan_from_spice():
    path_to_spice_folder = "./test_data/exp424"
    scan = Scan.from_spice(path_to_spice_folder, scan_num=34)
    assert scan.name == "scan0034"
    assert scan.scan_info.exp_id == "IPTS32124_CG4C_exp0424"
    assert scan.scan_info.scan_num == 34
    assert scan.scan_info.scan_title == "003 scan at Q=[0 0 2.5+0.8]"
    assert scan.sample_ub_info.sample_name == "NiTiO3"
    assert scan.sample_ub_info.type == "crystal"
    assert np.allclose(
        scan.sample_ub_info.ub_matrix,
        np.array(
            [[-0.016965, -0.026212, -0.071913], [-0.201388, -0.193307, 0.007769], [-0.108415, 0.1206, -0.003178]],
        ),
    )
    assert np.allclose(scan.sample_ub_info.u, [2.65493, 2.34763, -1.48203])
    assert np.allclose(scan.sample_ub_info.v, [-0.09782, -0.44945, -13.71879])
    assert scan.instrument_info.monochromator == "PG002"
    assert np.allclose(scan.instrument_info.collimation, (48, 40, 40, 120))
    assert len(scan.data) == 55
    assert len(scan.data["Pt."]) == 40
    assert np.allclose(scan.data["detector"][0:3], [569, 194, 40])


def test_scan_from_spice_hb1():
    path_to_spice_folder = "./test_data/IPTS31591_HB1_exp0917/exp917"
    Scan.from_spice(path_to_spice_folder, scan_num=45).plot()
    path_to_spice_folder = "./test_data/IPTS31591_HB1_exp0917/exp1111"
    for i in range(1, 5):
        Scan.from_spice(path_to_spice_folder, scan_num=i).plot()
    plt.show()


def test_scan_from_nexus():
    path_to_nexus_entry = "./test_data/IPTS32124_CG4C_exp0424/scan0034.h5"
    scan = Scan.from_nexus(path_to_nexus_entry)
    assert scan.name == "scan0034"
    assert scan.scan_info.scan_num == 34
    assert scan.scan_info.scan_title == "003 scan at Q=[0 0 2.5+0.8]"
    assert scan.sample_ub_info.sample_name == "NiTiO3"
    assert scan.sample_ub_info.type == "crystal"
    assert np.allclose(
        scan.sample_ub_info.ub_matrix,
        np.array(
            [[-0.016965, -0.026212, -0.071913], [-0.201388, -0.193307, 0.007769], [-0.108415, 0.1206, -0.003178]],
        ),
    )
    assert np.allclose(scan.sample_ub_info.u, [2.65493, 2.34763, -1.48203])
    assert np.allclose(scan.sample_ub_info.v, [-0.09782, -0.44945, -13.71879])
    assert scan.instrument_info.monochromator == "PG002"
    # assert np.allclose(scan.instrument_info.collimation, (48, 40, 40, 120))
    assert np.allclose(scan.instrument_info.collimation, (40, 100, 80, 120))
    assert len(scan.data) == 55
    assert len(scan.data["Pt."]) == 40
    assert np.allclose(scan.data["detector"][0:3], [569, 194, 40])


def test_get_data_columns():
    path_to_spice_folder = "./test_data/exp424"
    scan = Scan.from_spice(path_to_spice_folder, scan_num=42)
    assert len(scan.get_data_columns().keys()) == 55


def test_get_scan_data():
    path_to_spice_folder = "./test_data/exp424"
    scan = Scan.from_spice(path_to_spice_folder, scan_num=42)

    scan_data_1d = scan.get_data()

    x_data = np.arange(0.1, 4.1, 0.1)
    y_data = np.array(
        [481, 170, 22, 6, 2, 2, 1, 1, 2, 0, 0, 1, 4, 0, 2, 2, 6, 1, 1, 1]
        + [3, 1, 2, 1, 0, 2, 2, 1, 1, 2, 7, 14, 32, 76, 88, 101, 56, 34, 6, 5]
    )
    yerr_data = np.sqrt(y_data)
    mcu_data = np.array([60.0] * 40)
    time_data = np.array(
        [62.624, 66.241, 67.689, 65.918, 64.531, 65.345, 65.312, 65.053, 63.564, 62.493]
        + [61.6, 61.761, 62.119, 63.473, 64.994, 65.828, 66.575, 66.734, 67.131, 67.755]
        + [68.89, 69.928, 70.856, 71.629, 72.5, 72.871, 73.419, 73.649, 73.453, 73.78]
        + [73.738, 73.991, 74.261, 74.6, 74.638, 74.627, 75.343, 75.293, 75.183, 75.37]
    )

    assert np.allclose(scan_data_1d.x, x_data)
    assert np.allclose(scan_data_1d.y, y_data)
    assert np.allclose(scan_data_1d.err, yerr_data)
    assert scan.scan_info.preset_channel == "mcu"
    assert np.allclose(scan.data["mcu"], mcu_data)
    assert np.allclose(scan.data["time"], time_data)


def test_get_scan_data_norm():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    scan = Scan.from_nexus(nexus_file_name)
    scan_data_1d = scan.get_data(norm_to=(5, "mcu"))

    x_data = np.arange(0.1, 4.1, 0.1)
    y_data = np.array(
        [481, 170, 22, 6, 2, 2, 1, 1, 2, 0, 0, 1, 4, 0, 2, 2, 6, 1, 1, 1]
        + [3, 1, 2, 1, 0, 2, 2, 1, 1, 2, 7, 14, 32, 76, 88, 101, 56, 34, 6, 5]
    )
    yerr_data = np.sqrt(y_data)

    assert np.allclose(scan_data_1d.x, x_data)
    assert np.allclose(scan_data_1d.y, y_data / 12)
    assert np.allclose(scan_data_1d.err, yerr_data / 12)


def test_plot_scan_from_nexus():
    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    s1 = Scan.from_nexus(nexus_file_name)
    scan_data_1d = s1.get_data(norm_to=(30, "mcu"))
    assert scan_data_1d.label == "#42 "
    s1.plot(norm_to=(30, "mcu"))
    plt.show()
