import numpy as np

from tavi.data.spice_reader import read_spice_datafile, read_spice_ubconf


def test_read_spice_datafile_regular():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    (metadata, data, unrecognized, error_msg) = read_spice_datafile(path_to_spice_data)
    assert len(metadata) == 32
    assert metadata["scan"] == "34"
    assert data["Pt."].size == 40
    assert len(data) == 55


def test_read_spice_datafile_single_point():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0062.dat"
    (_, data, *_) = read_spice_datafile(path_to_spice_data)
    assert len(data) == 55
    assert data["Pt."].size == 1


def test_read_spice_datafile_no_ending():
    path_to_spice_data = "./test_data/exp416/Datafiles/CG4C_exp0416_scan0050.dat"
    (metadata, data, unrecognized, error_msg) = read_spice_datafile(path_to_spice_data)
    assert len(metadata) == 28


def test_read_spice_datafile_with_motor_errors():
    path_to_spice_data = "./test_data/exp815/Datafiles/HB1_exp0815_scan0001.dat"
    (*_, error_messages) = read_spice_datafile(path_to_spice_data)
    assert len(error_messages) == 44


def test_hb1_countfile_parsing():
    path_to_spice_data = "./test_data/IPTS31591_HB1_exp0917/exp1111/Datafiles/HB1_exp1111_scan0002.dat"
    (metadata, data, unrecognized, error_msg) = read_spice_datafile(path_to_spice_data)
    assert metadata["norm_channels"] == ["mcu_1", "mcu_2", "mcu_3"]
    assert metadata["labels"] == ["Px SF", "Px NSF", "Py SF"]


def test_read_spice_ubconf():
    path_to_spice_data = "./test_data/exp424/UBConf/UB02Jul2024_14108PM.ini"
    ubconf = read_spice_ubconf(path_to_spice_data)
    assert ubconf["UBMode"] == 1
    assert ubconf["AngleMode"] == 0
    assert ubconf["Energy"] == 4.8
    assert ubconf["UpperArc"] == 0
    assert np.allclose(
        ubconf["UBMatrix"],
        [-0.198255, -0.198255, 0.000000, 0.000000, 0.000000, -0.072359, 0.114463, -0.114463, -0.000000],
    )


def test_read_spice_ubconf_4_circle_html():
    path_to_spice_data = "./test_data/IPTS33347_HB1A_exp1046/exp1046/UBConf/UB20Feb2025_1645PM.ini"
    ubconf = read_spice_ubconf(path_to_spice_data)
    assert np.allclose(
        ubconf["UBMatrix"],
        [
            -0.00240830734807598,
            -0.0148785420823486,
            0.0404848827465389,
            0.201921571073141,
            0.192668717780335,
            0.00168620753782269,
            -0.106291529932909,
            0.121381235447627,
            0.00228599244682446,
        ],
    )
