"""Scan data model."""

from pydantic import BaseModel

from tavi.library.data.scan.metadata import ScanMetadata
from tavi.library.data.scan.ub import ScanUb
from tavi.library.data.scan.values import ScanValues


class Scan(BaseModel):
    """Scan data model."""

    ipts: str
    meta: ScanMetadata
    values: ScanValues
    ub: ScanUb
