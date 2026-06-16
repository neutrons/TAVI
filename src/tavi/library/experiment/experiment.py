"""
Experiment class.

Handles experimental data intake, extracting peak center, width etc.
"""

from pathlib import Path
from typing import Any

import numpy as np

from tavi.library.data.scan import RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.experiment.enum import FixedEnergyMode
from tavi.library.experiment.peak import DataPoint
from tavi.library.experiment.utilities import spice_to_mantid
from tavi.library.fit.fit import FitPackage, ModelName
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader
from tavi.library.storage.local_file_store import LocalFileStore


class Experiment:
    """Experiment class."""

    def __init__(
        self,
        mode: FixedEnergyMode = FixedEnergyMode.FIX_Ef,
        ei_or_ef: float = 0,
        loader: AbstractLoader = ORNLSpiceLoader(LocalFileStore()),
    ) -> None:
        """Init."""
        self.tavi_data: TaviData = TaviData(raw_scans={})
        self.mode = mode
        self.loader = loader
        self.set_ei_or_ef(ei_or_ef)

    def load_file(self, file_path: str) -> None:
        """
        Load a single scan.

        Directly using ORNL loader for now. will implement rounded approach later.
        """
        scan = self.loader.load(file_path)
        self.tavi_data.raw_scans[scan.uuid] = scan

    def load_folder(self, folder_path: str) -> None:
        """Load all SPICE .dat files in a folder of scans."""
        for file_path in sorted(Path(folder_path).glob("*.dat")):
            self.load_file(str(file_path))

    def get_hkl(
        self, scan_identifier: dict, use_title: bool = True, model_dict: list[tuple[ModelName, dict[str, Any]]] = []
    ) -> np.ndarray:
        """
        Extract the (h, k, l) from a scan title, rounded to 2 decimals.

        e.g. "scan_title =  (1.000019 -0.000008 0.499983) th4th, T = 4.3406 K"
             -> array([ 1.  , -0.  ,  0.5 ])
        """
        match self.loader:
            case ORNLSpiceLoader():
                scan_num = scan_identifier["scan_num"]
                IPTS = scan_identifier.get("IPTS", None)
                exp_num = scan_identifier.get("exp_num", None)
                return self.loader.get_hkl(
                    self.tavi_data, scan_num, IPTS, exp_num, use_title=use_title, model_dict=model_dict
                )
            case _:
                raise ValueError("Loader not implemented.")

    def get_peak_center(
        self, scan_identifier: dict, fit_package: FitPackage, model_dict: list[tuple[ModelName, dict[str, Any]]]
    ) -> DataPoint:
        """Find the center of the peak. It's used in refining UB matrix, which can be compared with SPICE results for validation."""
        match self.loader:
            case ORNLSpiceLoader():
                scan_num = scan_identifier["scan_num"]
                IPTS = scan_identifier.get("IPTS", None)
                exp_num = scan_identifier.get("exp_num", None)
                if self.mode is FixedEnergyMode.FIX_Ef:
                    ei_or_ef = self.ef
                else:
                    ei_or_ef = self.ei
                return self.loader.get_peak_center(
                    tavi_data=self.tavi_data,
                    scan_num=scan_num,
                    IPTS=IPTS,
                    exp_num=exp_num,
                    mode=self.mode,
                    ei_or_ef=ei_or_ef,
                    fit_package=fit_package,
                    model_dict=model_dict,
                )
            case _:
                raise ValueError("Loader not implemented.")

    def get_delta_q(self, scan_identifier: dict) -> np.ndarray:
        """Get delta q of a scan."""
        match self.loader:
            case ORNLSpiceLoader():
                scan_num = scan_identifier["scan_num"]
                IPTS = scan_identifier.get("IPTS", None)
                exp_num = scan_identifier.get("exp_num", None)
                if self.mode is FixedEnergyMode.FIX_Ef:
                    ei_or_ef = self.ef
                else:
                    ei_or_ef = self.ei
                return self.loader.get_delta_q(self.tavi_data, scan_num, IPTS, exp_num, self.mode, ei_or_ef)
            case _:
                raise ValueError("Loader not implemented.")

    def get_data_from_scan_number(self, scan_identifier: dict) -> RawScan:
        """
        Get scan object from a scan number.

        Args:
            scan_identifier: Dict with "scan_num" and optional "IPTS" / "exp_num"
                keys used to locate the scan.

        Returns:
            The matching RawScan.

        Raises:
            ValueError: If zero or more than one scan matches.

        """
        match self.loader:
            case ORNLSpiceLoader():
                scan_num = scan_identifier["scan_num"]
                IPTS = scan_identifier.get("IPTS", None)
                exp_num = scan_identifier.get("exp_num", None)
                return self.loader.get_data_from_scan_number(self.tavi_data, scan_num, IPTS, exp_num)
            case _:
                raise ValueError("Loader not implemented.")

    def get_two_theta(self, q_norm: float, ei: float, ef: float) -> float:
        """Get two_theta, only q_norm is required."""
        match self.loader:
            case ORNLSpiceLoader():
                return self.loader.get_two_theta(q_norm, ei, ef)
            case _:
                raise ValueError("Loader not implemented.")

    def get_psi(self, q_norm: float, ei: float, ef: float) -> float:
        """Get psi. Angle between ki and Q."""
        match self.loader:
            case ORNLSpiceLoader():
                return self.loader.get_psi(q_norm, ei, ef)
            case _:
                raise ValueError("Loader not implemented.")

    def set_ei_or_ef(self, e: float) -> None:
        """Set ei or ef based on mode."""
        if self.mode is FixedEnergyMode.FIX_Ef:
            self.ef = e
        else:
            self.ei = e

    def get_ei_ef(self, e: float) -> tuple[float, float]:
        """Get (ei, ef) given the complementary energy."""
        match self.loader:
            case ORNLSpiceLoader():
                if self.mode is FixedEnergyMode.FIX_Ef:
                    return self.loader.get_ei_ef(e, self.mode, self.ef)
                else:
                    return self.loader.get_ei_ef(e, self.mode, self.ei)
            case _:
                raise ValueError("Loader not implemented.")

    def create_sample(self, ub_path: str) -> Sample:
        """Can create a sample from a specific ub file if exist."""
        match self.loader:
            case ORNLSpiceLoader():
                ub = self.loader.load_hb1a_4c_ub(ub_path)
                sample = Sample(
                    OrientedLattice(
                        a=ub["a"], b=ub["b"], c=ub["c"], alpha=ub["alpha"], beta=ub["beta"], gamma=ub["gamma"]
                    )
                )
                sample.ol.UB = spice_to_mantid(ub["ub"])
                return sample
            case _:
                raise ValueError("loader not implemented.")
