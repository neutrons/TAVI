"""
Experiment class.

Handles experimental data intake, extracting peak center, width etc.
"""

from pathlib import Path

import numpy as np

from tavi.library.data.scan import RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.experiment.enum import FixedEnergyMode
from tavi.library.experiment.peak import DataPoint
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader
from tavi.library.storage.local_file_store import LocalFileStore


class Experiment:
    """Experiment class."""

    def __init__(
        self,
        mode: FixedEnergyMode = FixedEnergyMode.FIX_Ef,
        ei_or_ef: float = 0,
        loader: ORNLSpiceLoader = ORNLSpiceLoader(LocalFileStore()),
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

    def get_hkl(self, scan_identifier: dict) -> np.ndarray:
        """
        Extract the (h, k, l) from a scan title, rounded to 2 decimals.

        e.g. "scan_title =  (1.000019 -0.000008 0.499983) th4th, T = 4.3406 K"
             -> array([ 1.  , -0.  ,  0.5 ])
        """
        if isinstance(self.loader, ORNLSpiceLoader):
            scan_num = scan_identifier["scan_num"]
            IPTS = scan_identifier.get("IPTS", None)
            exp_num = scan_identifier.get("exp_num", None)
            return self.loader.get_hkl_from_title(self.tavi_data, scan_num, IPTS, exp_num)
        else:
            raise ValueError("Loader not implemented.")

    def create_peaks(self, scan_identifier: dict) -> DataPoint:
        """Create peak or peaks from scan numbers."""
        if isinstance(self.loader, ORNLSpiceLoader):
            scan_num = scan_identifier["scan_num"]
            IPTS = scan_identifier.get("IPTS", None)
            exp_num = scan_identifier.get("exp_num", None)
            if self.mode is FixedEnergyMode.FIX_Ef:
                ei_or_ef = self.ef
            else:
                ei_or_ef = self.ei
            return self.loader.create_peaks(self.tavi_data, scan_num, IPTS, exp_num, self.mode, ei_or_ef)
        else:
            raise ValueError("Loader not implemented.")

    def get_motor_angles(self, scan_identifier: dict) -> np.ndarray:
        """Generate Motor position from scan."""
        if isinstance(self.loader, ORNLSpiceLoader):
            scan_num = scan_identifier["scan_num"]
            IPTS = scan_identifier.get("IPTS", None)
            exp_num = scan_identifier.get("exp_num", None)
            return self.loader.get_motor_angles(self.tavi_data, scan_num, IPTS, exp_num)
        else:
            raise ValueError("Loader not implemented.")

    def get_delta_q(self, scan_identifier: dict) -> np.ndarray:
        """Get delta q of a scan."""
        if isinstance(self.loader, ORNLSpiceLoader):
            scan_num = scan_identifier["scan_num"]
            IPTS = scan_identifier.get("IPTS", None)
            exp_num = scan_identifier.get("exp_num", None)
            if self.mode is FixedEnergyMode.FIX_Ef:
                ei_or_ef = self.ef
            else:
                ei_or_ef = self.ei
            return self.loader.get_delta_q(self.tavi_data, scan_num, IPTS, exp_num, self.mode, ei_or_ef)
        else:
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
        if isinstance(self.loader, ORNLSpiceLoader):
            scan_num = scan_identifier["scan_num"]
            IPTS = scan_identifier.get("IPTS", None)
            exp_num = scan_identifier.get("exp_num", None)
            return self.loader.get_data_from_scan_number(self.tavi_data, scan_num, IPTS, exp_num)
        else:
            raise ValueError("Loader not implemented.")

    def get_two_theta(self, q_norm: float, ei: float, ef: float) -> float:
        """Get two_theta, only q_norm is required."""
        if isinstance(self.loader, ORNLSpiceLoader):
            return self.loader.get_two_theta(q_norm, ei, ef)

    def get_psi(self, q_norm: float, ei: float, ef: float) -> float:
        """Get psi. Angle between ki and Q."""
        if isinstance(self.loader, ORNLSpiceLoader):
            return self.loader.get_psi(q_norm, ei, ef)

    def set_ei_or_ef(self, e: float) -> None:
        """Set ei or ef based on mode."""
        if self.mode is FixedEnergyMode.FIX_Ef:
            self.ef = e
        else:
            self.ei = e

    def get_ei_ef(self, e: float) -> tuple[float, float]:
        """Get (ei, ef) given the complementary energy."""
        if isinstance(self.loader, ORNLSpiceLoader):
            if self.mode is FixedEnergyMode.FIX_Ef:
                return self.loader.get_ei_ef(e, self.mode, self.ef)
            else:
                return self.loader.get_ei_ef(e, self.mode, self.ei)
