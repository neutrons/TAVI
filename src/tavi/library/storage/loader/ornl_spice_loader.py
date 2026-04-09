"""ORNL Spice format loader."""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np

from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import UUID, Provenance, RawScan, Scan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.storage.interface.file_store_interface import FileStoreInterface
from tavi.library.storage.loader.interface.base import AbstractLoader


class ORNLSpiceLoader(AbstractLoader):
    """Loader for ORNL Spice format scan files."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """Initialize ORNL Spice loader with classifier."""
        super().__init__(filestore)
        self.classifier = RuleBasedClassifier(filestore)
        self.classification_rules = ORNLSpiceRuleSet()

    def load(self, file_path: str) -> Scan:
        """Load scan data."""
        uuid = UUID(file_path)
        values = self.parse_scan_values(file_path)
        meta = self.parse_metadata(file_path)
        tavi_meta = self.parse_tavi_metadata(file_path)
        prov = self.create_provenance(file_path)
        return self.adapt_scan_data(uuid=uuid, values=values, meta=meta, tavi_meta=tavi_meta, prov=prov)

    def get_scan_type(self) -> RawScanType:
        """Get scan type (ORNLSpice)."""
        return RawScanType.ORNLSpice

    def get_score(self, file_path: str) -> float:
        """Get score for scan."""
        return self.classifier.get_score(path, self.classification_rules)

    def parse_metadata(self, file_path: str) -> ScanMetadata:
        """Parse metadata."""
        with open(file_path, encoding="utf-8") as f:
            all_content = f.readlines()
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
        data = metadata | {"errors": error_messages} | {"others": others}
        return ScanMetadata(data)

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

        with open(file_path, encoding="utf-8") as f:
            all_content = f.readlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        metadata_list = headers[:index_col_name]

        preset_channel = ""
        preset_value = 0.0
        def_x, def_y = "", ""
        friendly_path = "IPTS"
        exp = "exp"
        s = "scan"
        for metadata_entry in metadata_list:
            if metadata_entry.startswith("# scan"):
                _, val = metadata_entry.split("=")
                scan += val.zfill(4)
            if metadata_entry.startswith("# proposal"):
                _, val = metadata_entry.split("=")
                friendly_path += val
            if metadata_entry.startswith("# experiment_number"):
                _, val = metadata_entry.split("=")
                exp += val.zfill(4)
            if metadata_entry.startswith("# proposal"):
                _, val = metadata_entry.split("=")
                friendly_path = val
            if metadata_entry.startswith("# preset_channel"):
                _, val = metadata_entry.split("=")
                preset_channel = val
            if metadata_entry.startswith("# preset_value"):
                _, val = metadata_entry.split("=")
                preset_value = val
            if metadata_entry.startswith("# def_x"):
                _, val = metadata_entry.split("=")
                default_x = val
            if metadata_entry.startswith("# def_y"):
                _, val = metadata_entry.split("=")
                def_y = val
        friendly_name = instrument_name + "_" + exp + "_" + s
        return TaviMetadata(
            default_axis=(def_x, def_y),
            friendly_name=friendly_name,
            friendly_path=friendly_path,
            normalization=(preset_channel, preset_value),
        )

    def parse_scan_values(self, file_path: str) -> ScanData:
        """Parse scan values."""
        with open(file_path, encoding="utf-8") as f:
            all_content = f.readlines()
        headers = [line.strip() for line in all_content if "#" in line]
        index_col_name = headers.index("# col_headers =")
        col_names = headers[index_col_name + 1].strip("#").split()
        col_values = np.genfromtxt(file_path, comments="#")
        data = dict()
        for col_name in col_names:
            # guard against invalid format
            if col_name[0].isdigit():  # can't start with digit, replace with _
                col_name = "_" + col_name
            attr_name = (
                col_name.replace("-", "_").replace(" ", "_").replace(".", "")
            )  # replace "-", " ", with "_", remove any "."
            try:
                data[attr_name] = col_values[:, col_names.index(col_name)]
            # sometimes data only have 1 entry, then we don't need to slice the data.
            except IndexError:
                data[attr_name] = np.array([col_values[col_names.index(col_name)]])
        return ScanData(data)

    def parse_external_metadata(self, file_path: str, ub_name: str) -> dict[str, Any]:
        """Parse corresponding file in ubconf as external metadata."""
        ubconf_path = Path(file_path).parent.parent.joinpath("UBConf").joinpath(ub_name)
        ubconf: dict[str, Any] = {}
        try:
            with open(ubconf_path, "r", encoding="utf-8") as f:
                all_content = f.readlines()
            if all_content[0] == "[UBMode]\n":
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
                    ub_matrix = matrix.attrib["matrix"].split(" ")
                ubconf.update({"UBMatrix": np.array([float(ub_matrix[i]) for i in range(9)])})
            return ubconf
        except FileNotFoundError:
            return {}

    def adapt_scan_data(
        self, uuid: UUID, values: ScanData, meta: ScanMetadata, tavi_meta: TaviMetadata, prov: Provenance
    ) -> RawScan:
        """Adapt scan data."""
        return RawScan(uuid=uuid, data=values, metadata=meta, tavimeta=tavi_meta, prov=prov)

    def create_provenance(self, file_path: str) -> Provenance:
        """Create provenance of the scan file."""
        uuid = UUID(file_path)
        weight = 1
        raw_file = file_path
        return Provenance(raw_file=raw_file, contributing_scans={uuid, weight})


if __name__ == "__main__":
    from pathlib import Path

    path = Path("/home/qmc/Neutrons/TAVI/test_data/exp815/Datafile/HB1_exp0815_scan001.dat")
    print(path.parent.parent.joinpath("UBConf"))
