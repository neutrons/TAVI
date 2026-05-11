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

    model_config = ConfigDict(arbitrary_types_allowed=True)
