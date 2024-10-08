import numpy as np

from tavi.data.scan import Scan


def test_scan_from_spice():
    path_to_spice_folder = "./test_data/exp424"
    scan = Scan.from_spice(path_to_spice_folder, scan_num=34)
    assert scan.name == "scan0034"
    assert scan.scan_info.scan_num == 34
    assert scan.scan_info.scan_title == "003 scan at Q=[0 0 2.5+0.8]"
    assert scan.sample_ub_info.sample_name == "NiTiO3"
    assert scan.sample_ub_info.type == "crystal"
    assert np.allclose(scan.sample_ub_info.u, [2.65493, 2.34763, -1.48203])
    assert np.allclose(scan.sample_ub_info.v, [-0.09782, -0.44945, -13.71879])
    assert scan.instrument_info.monochromator == "PG002"
    assert np.allclose(scan.instrument_info.collimation, (48, 40, 40, 120))
    assert len(scan.data) == 55
    assert len(scan.data["Pt."]) == 40
    assert np.allclose(scan.data["detector"][0:3], [569, 194, 40])


def test_scan_from_nexus():
    path_to_nexus_entry = "./test_data/spice_to_nxdict_test_scan0034.h5"
    scan = Scan.from_nexus(path_to_nexus_entry)
    assert scan.name == "scan0034"
    assert scan.name == "scan0034"
    assert scan.scan_info.scan_num == 34
    assert scan.scan_info.scan_title == "003 scan at Q=[0 0 2.5+0.8]"
    assert scan.sample_ub_info.sample_name == "NiTiO3"
    assert scan.sample_ub_info.type == "crystal"
    assert np.allclose(scan.sample_ub_info.u, [2.65493, 2.34763, -1.48203])
    assert np.allclose(scan.sample_ub_info.v, [-0.09782, -0.44945, -13.71879])
    assert scan.instrument_info.monochromator == "PG002"
    assert np.allclose(scan.instrument_info.collimation, (48, 40, 40, 120))
    assert len(scan.data) == 55
    assert len(scan.data["Pt."]) == 40
    assert np.allclose(scan.data["detector"][0:3], [569, 194, 40])
