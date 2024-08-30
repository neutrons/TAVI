# -*- coding: utf-8 -*-

from tavi.data.scan import Scan


def test_load_scan():

    nexus_file_name = "./test_data/IPTS32124_CG4C_exp0424/scan0042.h5"
    scan = Scan.from_nexus(nexus_file_name)
