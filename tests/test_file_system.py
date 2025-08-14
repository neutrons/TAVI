import os
from dataclasses import fields

import numpy as np

from tavi.tavi_model.FileSystem.dataloader import load_folder


def test_load_raw_data():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    scan = tavi_project.scans["CG4C_exp0424_scan0073.dat"]
    assert len(scan.data.Pt) == len(scan.data.q)
    assert len(scan.data.Pt) == len(scan.data.h)
    assert len(scan.data.Pt) == len(scan.data.detector)
    assert scan.data.detector[0:3].all() == np.array([2.608e03, 4.360e02, 4.600e01]).all()
    assert scan.data.focal_length[-1] == 72.6996
    assert len(tavi_project.scans["CG4C_exp0424_scan0073.dat"].error_message) == 0


def test_load_all_raw_data_are_pased_correctly():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    scan = tavi_project.scans["CG4C_exp0424_scan0001.dat"]
    for field in fields(scan.data):
        assert getattr(scan.data, field.name) is not None


def test_load_all_raw_meata_data_are_pased_correctly():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    scan = tavi_project.scans["CG4C_exp0424_scan0041.dat"]
    for field in fields(scan.metadata):
        assert getattr(scan.metadata, field.name) is not None


def test_load_raw_data_with_error():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    scan = tavi_project.scans["CG4C_exp0424_scan0041.dat"]
    assert len(scan.data.Pt) == 1
    assert len(scan.error_message) == 4
    assert scan.metadata.Sum_of_Counts == "0"


def test_all_ub_files_are_loaded():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    ub = tavi_project.ubconf
    assert len(ub) == 13


def test_load_ub():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    ub = tavi_project.ubconf["UB02Jul2024_14108PM.ini"]
    assert ub.UBMode == 1
    assert ub.Peak1.all() == np.array([[0.0, 0.0, 3.0, 53.24, 32.865, 2.3108, 0.0, 4.799999, 4.799998]]).all()


def test_load_ub_tmp():
    folder_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(folder_dir, os.pardir, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    ub_tmp = tavi_project.ubconf["tmp-UB02Jul2024_21029PM.ini"]
    assert ub_tmp.UBMode == 1
    assert ub_tmp.LatticeParams.all() == np.array([5.044, 5.044, 13.82, 90.0, 90.0, 120.0]).all()
