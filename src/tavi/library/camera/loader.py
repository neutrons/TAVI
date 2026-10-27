"""Loader for DEMAND/HB3A Anger camera SPICE XML point files and their UB matrices."""

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from tavi.library.experiment.utilities import spice_to_mantid

_SHAPE = re.compile(r"^(\w+?)(?:\[([\d,]+)\])?$")
_CASTS = {"FLOAT32": float, "FLOAT64": float, "INT32": int, "INT64": int}


class AngerCameraLoader:
    """Parse Anger camera SPICE XML points and UB matrices into one in-memory store."""

    def __init__(self) -> None:
        """Start with an empty store; parsing accumulates points keyed by filename."""
        self.data: dict = {"ub": {}}

    def _value(self, elem: ET.Element) -> str | int | float | np.ndarray:
        """Convert one SPICE XML leaf to a scalar, array, or image per its `type` attribute."""
        text = (elem.text or "").strip()
        match = _SHAPE.match(elem.get("type", ""))
        cast = _CASTS.get(match.group(1)) if match else None
        if cast is None:  # untyped: Header strings, DetectorDataFilePath
            return text
        if match.group(2) is None:
            return cast(text)
        shape = tuple(int(n) for n in match.group(2).split(","))
        return np.fromstring(text, dtype=np.int64 if cast is int else np.float64, sep="\t").reshape(shape)

    def _walk(self, elem: ET.Element) -> dict | str | int | float | np.ndarray:
        """Recurse into section elements; leaves become typed values."""
        if len(elem) == 0:
            return self._value(elem)
        return {child.tag: self._walk(child) for child in elem}

    def parse_spice_xml(self, path: str | Path) -> None:
        """Parse one DEMAND/HB3A SPICE XML point file into self.data under its filename."""
        path = Path(path)
        root = ET.parse(path).getroot()
        point = {"file": dict(root.attrib)} | {child.tag: self._walk(child) for child in root}
        self.data[path.name] = point

    def parse_folder_xml(self, path: str | Path) -> None:
        """Parse every XML file in a folder into self.data, skipping anything else."""
        for entry in sorted(Path(path).iterdir()):
            if entry.is_file() and entry.suffix.lower() == ".xml":
                self.parse_spice_xml(entry)

    def load_ub(self, path: str | Path) -> None:
        """Load the 3x3 UB matrix from a SPICE UBConf .dat file into self.ub, in SPICE convention."""
        path = Path(path)
        root = ET.parse(path).getroot()
        self.data["ub"][path.name] = spice_to_mantid(
            np.fromstring(root.find("matrix").get("matrix"), dtype=np.float64, sep=" ").reshape(3, 3), version="old"
        )

    def set_ub(self, ub: str, ub_name: str) -> None:
        """Store a UB matrix given as nine whitespace-separated numbers in SPICE convention."""
        self.data["ub"][ub_name] = spice_to_mantid(
            np.fromstring(ub, dtype=np.float64, sep=" ").reshape(3, 3), version="old"
        )
