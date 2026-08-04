"""ORNL Spice format loader."""

import logging
import re
import warnings
import xml.etree.ElementTree as ET
from typing import Any, Optional

import numpy as np

from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import UUID, Provenance, RawScan, Scan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.data.tavi_data import TaviData
from tavi.library.experiment.enum import FixedEnergyMode
from tavi.library.experiment.peak import DataPoint, MotorAngles
from tavi.library.experiment.utilities import SE2K, get_angle_from_triangle, get_side_from_triangle
from tavi.library.fit import Fit, FitPackage, ModelName
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.interface.base import AbstractLoader

logger = logging.getLogger(__name__)


class ORNLSpiceLoader(AbstractLoader):
    """Loader for ORNL Spice format scan files."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize ORNL Spice loader with classifier."""
        super().__init__(filestore)
        self.classifier = RuleBasedClassifier(filestore)
        self.classification_rules = ORNLSpiceRuleSet()

    def load(self, file_path: str) -> Scan:
        """Load scan data."""
        uuid = self.generate_uuid(file_path)
        values = self.parse_scan_values(file_path)
        meta = self.parse_metadata(file_path)
        tavi_meta = self.parse_tavi_metadata(file_path)
        prov = self.create_provenance(file_path)

        # get ubconf file name
        ub_name = meta.ubconf
        ubconf = self.parse_external_metadata(file_path, ub_name)
        # add it to MetaData's data entry
        meta.data["ORNL Metadata"].update(ubconf)
        return self.adapt_scan_data(uuid=uuid, values=values, meta=meta, tavi_meta=tavi_meta, prov=prov)

    def get_scan_type(self) -> RawScanType:
        """Get scan type (ORNLSpice)."""
        return RawScanType.ORNLSpice

    def get_score(self, file_path: str) -> float:
        """Get score for scan."""
        self.classifier.set_filestore(self.filestore)
        return self.classifier.get_score(file_path, self.classification_rules)

    def parse_metadata(self, file_path: str) -> ScanMetadata:
        """Parse metadata."""
        f = self.filestore.read_text_file(file_path=file_path)
        all_content = f.splitlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        col_names = headers[index_col_name + 1].strip("#").split()
        # remove the dot before it causes problem
        # index_of_pt = col_names.index("Pt.")
        # col_names[index_of_pt] = "Pt"
        metadata_list = headers[:index_col_name]
        error_messages = headers[index_col_name + 2 :]

        index_sum_count = [i for i, header in enumerate(headers) if header.startswith("# Sum of Counts =")]
        # in case "Sum of Counts" doesn't exist
        # happens to the last scan after beam is down
        if len(index_sum_count) != 0:
            metadata_list += headers[index_sum_count[0] :]
            error_messages = error_messages[: index_sum_count[0] - len(headers)]

        metadata = {}
        others = []

        for metadata_entry in metadata_list:
            line = metadata_entry.strip("# ")

            if "completed" in line or "stopped" in line:  # last line
                parts = line.split(" ")
                end_time = parts[3] + " " + parts[0] + " " + parts[1]
                metadata.update({"end_time": end_time})
            # elif line[-1] == "=":  # empty line
            #     unused.append(line[:-2])  # remove  " ="
            elif "=" in line:  # useful line
                parts = line.split("=")
                key = parts[0].strip()
                val = "=".join(parts[1:])[1:]  # remove the first space character
                metadata.update({key: val})
            else:  # how did you get here?
                others.append(line)

        if metadata.get("preset_type") == "countfile":  # HB1 in polarization mode
            countfile = []
            for metadata_entry in metadata_list:
                if metadata_entry.startswith("# countfile"):
                    _, val = metadata_entry.split("=")
                    countfile.append(val.strip())
            metadata.update({"countfile": ", ".join(countfile)})

        metadata.update({"errors": error_messages, "others": others})
        data = {"ORNL Metadata": metadata}
        categories = {"ORNL Metadata": "ORNL Metadata"}

        return ScanMetadata(data=data, categories=categories)

    def parse_tavi_metadata(self, file_path: str) -> TaviMetadata:
        """Parse metadata."""
        instrument_name = ""
        if "HB1A" in file_path:
            instrument_name = "HB1A"
        if "CG4C" in file_path:
            instrument_name = "CG4C"
        if "HB1" in file_path:
            instrument_name = "HB1"
        if "HB3" in file_path:
            instrument_name = "HB3"

        f = self.filestore.read_text_file(file_path=file_path)
        all_content = f.splitlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        metadata_list = headers[:index_col_name]

        preset_channel = ""
        preset_value = 0.0
        def_x, def_y = "", ""
        friendly_path = "IPTS-"
        exp = "exp"
        s = "scan"
        for metadata_entry in metadata_list:
            if metadata_entry.startswith("# scan "):
                _, val = metadata_entry.split("=")
                s += val.strip().zfill(4)
            if metadata_entry.startswith("# proposal"):
                _, val = metadata_entry.split("=")
                friendly_path += val.strip()
            if metadata_entry.startswith("# experiment_number"):
                _, val = metadata_entry.split("=")
                exp += val.strip().zfill(4)
            if metadata_entry.startswith("# preset_channel"):
                _, val = metadata_entry.split("=")
                preset_channel = val
            if metadata_entry.startswith("# preset_value"):
                _, val = metadata_entry.split("=")
                preset_value = val
            if metadata_entry.startswith("# def_x"):
                _, val = metadata_entry.split("=")
                def_x = val
            if metadata_entry.startswith("# def_y"):
                _, val = metadata_entry.split("=")
                def_y = val
        friendly_name = instrument_name + "_" + exp + "_" + s
        return TaviMetadata(
            default_axis=(def_x.strip(), def_y.strip()),
            friendly_name=friendly_name,
            friendly_path=friendly_path,
            normalization=(preset_channel, preset_value),
        )

    def parse_scan_values(self, file_path: str) -> ScanData:
        """Parse scan values."""
        f = self.filestore.read_text_file(file_path=file_path)
        all_content = f.splitlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        col_names = headers[index_col_name + 1].strip("#").split()
        try:
            with warnings.catch_warnings():
                # Treat all warnings as exceptions within this block
                warnings.simplefilter("error")
                try:
                    col_values = np.genfromtxt(file_path, comments="#")
                except ValueError as e:
                    logger.error(e)
                    col_values = np.array(None)
        except Warning as e:
            # exception happens when there is no valid measurements but all warnings.
            # see HB1_exp0815_scan0001.dat file
            logger.error(e)
            col_values = np.array(None)
        data = dict()
        for col_index, col_name in enumerate(col_names):
            # guard against invalid format
            if col_name[0].isdigit():  # can't start with digit, replace with _
                col_name = "_" + col_name
            attr_name = (
                col_name.replace("-", "_").replace(" ", "_").replace(".", "")
            )  # replace "-", " ", with "_", remove any "."
            if col_values.ndim > 1:
                data[attr_name] = col_values[:, col_index]
            # sometimes data only have 1 entry, then we don't need to slice the data.
            elif col_values.ndim == 1:
                data[attr_name] = np.array([col_values[col_index]])
            else:
                data[attr_name] = []
        return ScanData(data=data)

    def parse_external_metadata(self, file_path: str, ub_name: str) -> dict[str, Any]:
        """Parse corresponding file in ubconf as external metadata."""
        # ub_name is recorded by SPICE as a full Windows path (e.g.
        # "C:\\SPICE\\User\\exp1046\\UBConf\\file.dat"); keep only the filename
        # so it resolves against the local UBConf folder.
        ub_name = ub_name.replace("\\", "/").rsplit("/", 1)[-1]
        root_path = file_path
        for _ in range(2):
            root_path = self.filestore.get_parent(root_path)
        ubconf_path = self.filestore.join_path(root_path, "UBConf")
        ubconf_path = self.filestore.join_path(ubconf_path, ub_name)
        try:
            return self._parse_ubconf(ubconf_path=ubconf_path)
        except FileNotFoundError:
            return {}

    def _parse_ubconf(self, ubconf_path: str) -> dict[str, Any]:
        """Parse a .ini file in ubconf folder for ORNL TAS data."""
        ubconf: dict[str, Any] = {}
        f = self.filestore.read_text_file(file_path=ubconf_path)
        all_content = f.splitlines()
        if all_content[0] == "[UBMode]":
            for idx, line in enumerate(all_content):
                if line.strip() == "":
                    continue  # skip if empty
                elif line.strip()[0] == "[":
                    continue  # skiplines like "[xx]"
                key, val = line.strip().split("=")

                if key == "Mode":
                    mode_name = all_content[idx - 1].strip()
                    if mode_name == "[UBMode]":
                        ubconf.update({"UBMode": int(val)})
                    elif mode_name == "[AngleMode]":
                        ubconf.update({"AngleMode": int(val)})
                elif "," in val:  # string of vector to array
                    ubconf.update({key: np.array([float(v) for v in val.strip('"').split(",")])})
                elif val == '""':  # no value
                    pass
                else:  # float
                    ubconf.update({key: float(val)})
        else:  # xml junk from C#
            tree = ET.parse(ubconf_path)
            root = tree.getroot()
            for matrix in root.findall("matrix"):
                ub_matrix = matrix.attrib["matrix"].split()
            ubconf.update({"UBMatrix": np.array([float(ub_matrix[i]) for i in range(9)])})
        return ubconf

    def adapt_scan_data(
        self, uuid: UUID, values: ScanData, meta: ScanMetadata, tavi_meta: TaviMetadata, prov: Provenance
    ) -> RawScan:
        """Adapt scan data."""
        return RawScan(uuid=uuid, data=values, metadata=meta, tavimeta=tavi_meta, prov=prov)

    def create_provenance(self, file_path: str) -> Provenance:
        """Create provenance of the scan file."""
        uuid = self.generate_uuid(file_path)
        weight = 1
        raw_file = file_path
        return Provenance(raw_file=raw_file, contributing_scans={uuid: weight})

    # ----------------------------------------------------------------------------------------------
    #                                      ORNL specific
    # ----------------------------------------------------------------------------------------------
    # make static
    def get_hkl(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
        use_title: bool = True,
        model_dict: list[tuple[ModelName, dict[str, Any]]] = [],
    ) -> np.ndarray:
        """
        Extract the (h, k, l), user can define if they want to just get hkl from a scan title, rounded to 2 decimals. or use a fitting function to fit.

        e.g. "scan_title =  (1.000019 -0.000008 0.499983) th4th, T = 4.3406 K"
             -> array([ 1.  , -0.  ,  0.5 ])
        """
        scan = self.get_data_from_scan_number(tavi_data=tavi_data, scan_num=scan_num, IPTS=IPTS, exp_num=exp_num)
        if use_title:
            scan_title = scan.metadata.scan_title

            m = re.search(r"\(([^)]+)\)", scan_title)
            if m is None:
                raise ValueError(f"No (h k l) found in scan_title: {scan_title!r}")
            hkl = np.array([float(v) for v in re.split(r"[,\s]+", m.group(1).strip())])
            return np.round(hkl, 2)
        else:
            tol = 1e-3
            # Using fitting to find hkl
            fit = Fit(package=FitPackage.lmfit)
            def_y = scan.metadata.def_y
            y = scan.data.data[def_y]
            if abs(max(scan.data.h) - min(scan.data.h)) > tol:
                h_center = self._fit_centers(fit, scan.data.h, y, model_dict)
                k_center = np.mean(scan.data.k)
                l_center = np.mean(scan.data.l)
            elif abs(max(scan.data.k) - min(scan.data.k)) > tol:
                h_center = np.mean(scan.data.h)
                k_center = self._fit_centers(fit, scan.data.k, y, model_dict)
                l_center = np.mean(scan.data.l)
            elif abs(max(scan.data.l) - min(scan.data.l)) > tol:
                h_center = np.mean(scan.data.h)
                k_center = np.mean(scan.data.k)
                l_center = self._fit_centers(fit, scan.data.l, y, model_dict)
            else:
                h_center = np.mean(scan.data.h)
                k_center = np.mean(scan.data.k)
                l_center = np.mean(scan.data.l)
            hkl = (h_center, k_center, l_center)
            return hkl

    def get_peak_center(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
        mode: FixedEnergyMode = FixedEnergyMode.FIX_Ef,
        fixed_energy: float = 0,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: list[tuple[ModelName, dict[str, Any]]] = [],
    ) -> list[DataPoint]:
        """Create one DataPoint per fitted peak center in the scan."""
        scan = self.get_data_from_scan_number(tavi_data, scan_num, IPTS, exp_num)
        hkl = self.get_hkl(tavi_data, scan_num, IPTS, exp_num)
        motor_angles_list = self.fit_motor_angles(tavi_data, scan_num, IPTS, exp_num, fit_package, model_dict)
        # Detect if this is an inelastic scan
        fit = Fit(package=fit_package)
        y = scan.metadata.def_y
        if abs(max(scan.data.e) - min(scan.data.e)) > 0.1:
            e_center = self._fit_centers(fit, scan.data.e, scan.data.data[y], model_dict)
            eis, efs = [], []
            for i in range(len(e_center)):
                ei, ef = self.get_ei_ef(e_center[i], mode, fixed_energy)
                eis.append(ei)
                efs.append(ef)
            # One DataPoint per peak, pairing each motor-angle peak with its energy peak.
            return [
                DataPoint(hkl=hkl, ei=ei, ef=ef, angles=angles) for angles, ei, ef in zip(motor_angles_list, eis, efs)
            ]
        else:
            e_center = np.mean(scan.data.e)
        ei, ef = self.get_ei_ef(e_center, mode, fixed_energy)
        return [DataPoint(hkl=hkl, ei=ei, ef=ef, angles=angles) for angles in motor_angles_list]

    def fit_motor_angles(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: list[tuple[ModelName, dict[str, Any]]] = [],
    ) -> list[MotorAngles]:
        """
        Generate fitted motor positions, one MotorAngles per fitted peak center.

        The scanned (default-x) motor is fit with model_dict; when it has a paired
        motor (omega for a 2theta scan, s1 for an s2 scan) that motor is fit too and
        centers are paired by sorted order. Motors that are not scanned use their
        mean over the scan. A list with one entry per fitted center is returned.
        """
        scan = self.get_data_from_scan_number(tavi_data=tavi_data, scan_num=scan_num, IPTS=IPTS, exp_num=exp_num)
        def_x, def_y = scan.metadata.def_x, scan.metadata.def_y
        if def_x[0].isdigit():
            def_x = "_" + def_x
        if def_y[0].isdigit():
            def_y = "_" + def_y
        x, y = scan.data.data[def_x], scan.data.data[def_y]
        fit = Fit(package=fit_package)

        # fitted: motor -> per-peak centers; means: motor -> single value shared by all peaks.
        fitted: dict[str, list[float]] = {}
        means: dict[str, float] = {}
        match def_x:
            case "_2theta":
                fitted["two_theta"] = self._fit_centers(fit, x, y, model_dict)
                fitted["omega"] = self._fit_centers(fit, scan.data.omega, y, model_dict)
                means["chi"] = np.mean(scan.data.chi)
                means["phi"] = np.mean(scan.data.phi)
            case "omega":
                means["two_theta"] = np.mean(scan.data._2theta)
                fitted["omega"] = self._fit_centers(fit, x, y, model_dict)
                means["chi"] = np.mean(scan.data.chi)
                means["phi"] = np.mean(scan.data.phi)
            case "chi":
                means["two_theta"] = np.mean(scan.data._2theta)
                means["omega"] = np.mean(scan.data.omega)
                fitted["chi"] = self._fit_centers(fit, x, y, model_dict)
                means["phi"] = np.mean(scan.data.phi)
            case "phi":
                means["two_theta"] = np.mean(scan.data._2theta)
                means["omega"] = np.mean(scan.data.omega)
                means["chi"] = np.mean(scan.data.chi)
                fitted["phi"] = self._fit_centers(fit, x, y, model_dict)
            case "s2":
                fitted["s2"] = self._fit_centers(fit, x, y, model_dict)
                fitted["s1"] = self._fit_centers(fit, scan.data.s1, y, model_dict)
                means["sgl"] = np.mean(scan.data.sgl)
                means["sgu"] = np.mean(scan.data.sgu)
            case "s1":
                means["s2"] = np.mean(scan.data.s2)
                fitted["s1"] = self._fit_centers(fit, x, y, model_dict)
                means["sgl"] = np.mean(scan.data.sgl)
                means["sgu"] = np.mean(scan.data.sgu)
            case "sgl":
                means["s2"] = np.mean(scan.data.s2)
                means["s1"] = np.mean(scan.data.s1)
                fitted["sgl"] = self._fit_centers(fit, x, y, model_dict)
                means["sgu"] = np.mean(scan.data.sgu)
            case "sgu":
                means["s2"] = np.mean(scan.data.s2)
                means["s1"] = np.mean(scan.data.s1)
                means["sgl"] = np.mean(scan.data.sgl)
                fitted["sgu"] = self._fit_centers(fit, x, y, model_dict)
            case _:
                raise ValueError("Not implemented yet.")

        # All fitted motors share the same peak count; guard against a mismatch.
        n_peaks = min(len(centers) for centers in fitted.values())
        all_motors = ("two_theta", "omega", "chi", "phi", "s2", "s1", "sgl", "sgu")
        angles_list = []
        for i in range(n_peaks):
            angles_dict = {motor: (fitted[motor][i] if motor in fitted else means.get(motor)) for motor in all_motors}
            angles_list.append(MotorAngles(angles_dict=angles_dict))
        return angles_list

    @staticmethod
    def _fit_centers(
        fit: Fit, x: np.ndarray, y: np.ndarray, model_dict: list[tuple[ModelName, dict[str, Any]]]
    ) -> list[float]:
        """
        Fit (x, y) with model_dict and return every fitted center, sorted ascending.

        Components without a "center" parameter (e.g. background models) are ignored.
        """
        result = fit.fit(x, y, model_dict)
        centers = [comp.values["center"] for comp in result.components.values() if "center" in comp.values]
        return sorted(centers)

    def get_data_point_closest_to_center(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: list[tuple[ModelName, dict[str, Any]]] = [],
        mode: FixedEnergyMode = FixedEnergyMode.FIX_Ef,
        fixed_energy: float = 0,
    ) -> list:
        """Create a DataPoint from the scan row whose column value is closest to center."""
        scan = self.get_data_from_scan_number(tavi_data, scan_num, IPTS, exp_num)
        # get default x, y data
        def_x = scan.metadata.def_x
        def_y = scan.metadata.def_y
        if def_x[0].isdigit():
            def_x = "_" + def_x
        if def_y[0].isdigit():
            def_y = "_" + def_y

        if def_x in ["s1", "s2", "omega"]:
            x = self.get_delta_q(tavi_data, scan_num, IPTS, exp_num, mode, fixed_energy)
        else:
            x = np.asarray(scan.data.data[def_x])
        y = np.asanyarray(scan.data.data[def_y])
        fit = Fit(fit_package)
        # there may be more than 1 centers depending on user chosen fit models
        centers = self._fit_centers(fit, x, y, model_dict)
        # find the index whose column value is closest to center
        idx_list = [np.argmin(np.abs(x - center)) for center in centers]
        # now get ORNL specfici information to create a data point
        data_point_list = []
        for idx in idx_list:
            h, k, l = scan.data.h[idx], scan.data.k[idx], scan.data.l[idx]
            # Issue here
            try:
                ei = scan.data.ei[idx]
                ef = ei - scan.data.e[idx]
            except AttributeError:
                ef = scan.data.ef[idx]
                ei = ef - scan.data.e[idx]
            if "s2" in scan.data.data:
                s1 = scan.data.s1[idx]
                s2 = scan.data.s2[idx]
                sgl = scan.data.sgl[idx]
                sgu = scan.data.sgu[idx]
                data_point_list.append(
                    DataPoint(
                        (h, k, l),
                        ei,
                        ef,
                        MotorAngles(angles_dict={"two_theta": s2, "omega": s1, "sgl": sgl, "sgu": sgu}),
                    )
                )
            elif "_2theta" in scan.data.data:
                two_theta = scan.data._2theta[idx]
                omega = scan.data.omega[idx]
                chi = scan.data.chi[idx]
                phi = scan.data.phi[idx]
                data_point_list.append(
                    DataPoint(
                        (h, k, l),
                        ei,
                        ef,
                        MotorAngles(angles_dict={"two_theta": two_theta, "omega": omega, "chi": chi, "phi": phi}),
                    )
                )
            else:
                raise ValueError("Data not of ORNL type.")
        return data_point_list

    def get_delta_q(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
        mode: FixedEnergyMode = FixedEnergyMode.FIX_Ef,
        fixed_energy: float = 0,
    ) -> np.ndarray:
        """Get delta q of a scan."""
        scan = self.get_data_from_scan_number(tavi_data=tavi_data, scan_num=scan_num, IPTS=IPTS, exp_num=exp_num)
        try:
            qs = scan.data.q
        except AttributeError:
            es = scan.data.e
            if mode == FixedEnergyMode.FIX_Ef:
                ef = fixed_energy
                eis = [e + ef for e in es]
                efs = [ef] * len(eis)
            else:
                ei = fixed_energy
                efs = [e - ei for e in es]
                eis = [ei] * len(efs)
            kis = np.array([SE2K(ei) for ei in eis])
            kfs = np.array([SE2K(ef) for ef in efs])
            two_thetas = scan.data.s2 if "s2" in dir(scan.data) else scan.data._2theta
            qs = np.array(
                [
                    get_side_from_triangle(ki, kf, np.radians(two_theta))
                    for ki, kf, two_theta in zip(kis, kfs, two_thetas)
                ]
            )

        q_diff = np.max(qs) - np.min(qs)
        mid_idx = (len(qs) - 1) // 2

        if q_diff > 1.1e-4:  # q changing, must be a th2th scan
            return qs - qs[mid_idx]
        else:  # q not changing, must be a s1 scan
            q_abs = np.mean(qs)

            if "s1" in dir(scan.data):  # using "s1" by default
                angles = np.asanyarray(scan.data.s1)
            elif "omega" in dir(scan.data):
                angles = np.asanyarray(scan.data.omega)
            else:
                raise AttributeError("No s1 or omega in data. Can't calculate delta q")
            return np.radians(angles - angles[mid_idx]) * q_abs

    def get_data_from_scan_number(
        self,
        tavi_data: TaviData,
        scan_num: int,
        IPTS: Optional[int] = None,
        exp_num: Optional[int] = None,
    ) -> RawScan:
        """
        Get scan object from a scan number.

        Args:
            tavi_data: Container of loaded raw scans to search.
            scan_num: Scan number to look up.
            IPTS: Optional IPTS proposal number to disambiguate across experiments.
            exp_num: Optional experiment number to disambiguate within an IPTS.

        Returns:
            The matching RawScan.

        Raises:
            ValueError: If zero or more than one scan matches.

        """
        scan_tag = "scan" + str(scan_num).zfill(4)
        exp_tag = "exp" + str(exp_num).zfill(4) if exp_num is not None else None
        ipts_tag = "IPTS-" + str(IPTS) if IPTS is not None else None

        matches = []
        for raw_scan in tavi_data.raw_scans.values():
            name = raw_scan.tavimeta.friendly_name
            path = raw_scan.tavimeta.friendly_path
            if scan_tag not in name:
                continue
            if exp_tag is not None and exp_tag not in name:
                continue
            if ipts_tag is not None and ipts_tag not in path:
                continue
            matches.append(raw_scan)

        if len(matches) == 0:
            raise ValueError(f"No scan found with scan_num={scan_num}, IPTS={IPTS}, exp_num={exp_num}.")
        if len(matches) > 1:
            raise ValueError(
                f"{len(matches)} scans match scan_num={scan_num}. Specify IPTS and/or exp_num to disambiguate."
            )
        return matches[0]

    def get_two_theta(self, q_norm: float, ei: float, ef: float) -> float:
        """Get two_theta, only q_norm is required."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        two_theta_rad = get_angle_from_triangle(ki, kf, q_norm)
        return two_theta_rad

    def get_psi(self, q_norm: float, ei: float, ef: float) -> float:
        """Get psi. Angle between ki and Q."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        psi_rad = get_angle_from_triangle(ki, q_norm, kf)
        return psi_rad

    def get_ei_ef(self, e: float, mode: FixedEnergyMode, fixed_energy: float) -> tuple[float, float]:
        """Get (ei, ef) given the complementary energy."""
        if mode is FixedEnergyMode.FIX_Ef:
            ef = fixed_energy
            ei = e + ef
            return ei, ef
        else:
            ei = fixed_energy
            ef = ei - e
            return ei, ef

    def load_hb1a_4c_ub(self, file_path: str) -> dict:
        """C sharp ub specific to ORNL/HB1A in 4 circle mode."""
        root = ET.parse(file_path).getroot()

        cell = root.find("unitcell").attrib
        a, b, c = float(cell["a"]), float(cell["b"]), float(cell["c"])
        alpha, beta, gamma = float(cell["alpha"]), float(cell["beta"]), float(cell["gamma"])
        wavelength = float(root.find("wavelength").attrib["lambda"])
        ub = np.array([float(x) for x in root.find("matrix").attrib["matrix"].split()]).reshape(3, 3)
        return dict(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, ub=ub, wavelength=wavelength)
