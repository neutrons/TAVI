"""Plot data model."""

from typing import Optional

import numpy as np
from pydantic import BaseModel, ConfigDict

from tavi.library.data.scan import UUID, UUIDFactory


class Plot(BaseModel):
    """Represents a single plottable dataset derived from a scan."""

    uuid: UUID = UUIDFactory()
    x: np.ndarray
    y: np.ndarray
    err: np.ndarray
    scan_name: str
    normalized_by: Optional[str]
    x_name: str
    y_name: str
    error_name: str
    source_scan_uuid: UUID
    """uuid of the immutable RawScan this plot derived from. Regenerate a new Plot from it instead of mutating this one."""

    model_config = ConfigDict(arbitrary_types_allowed=True)
