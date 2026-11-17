"""What a group fit is asked to do, and what it produced for each scan."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from tavi.library.data.fit_entry import FitCurve, FitData, FitResultSummary, FitSpec
from tavi.library.data.scan import UUID


class GroupFitRequest(BaseModel):
    """
    Everything ``GroupFitModel`` needs to fit a group: each scan's raw panel spec, and the columns they share.

    Shaped like the Group Fit window's own state - ``specs`` is what ``GroupFitWindow.fit_specs()``
    returns - so the view hands over raw text and the model does all parsing, as ``FitRequest`` does.
    """

    specs: dict[UUID, FitSpec]
    """Scan uuid -> that scan's spec, in fit order. The first scan's spec is the group model; with
    ``chain`` on, a later scan overrides it only where the user changed a parameter for that scan."""
    y_col: str
    x_col: Optional[str] = None
    """``None`` fits each scan against its own ``tavimeta.default_axis[0]``, as ``GroupFitParams.x_col``."""
    chain: bool = True
    """Start each scan after the first from the previous scan's fitted values - see ``GroupFitParams.chain``."""
    averages: list[str] = []
    """Columns reduced to a per-scan mean and spread instead of fit, e.g. sample temperature."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


class SuggestGroupParamsRequest(BaseModel):
    """
    Which scan of a group to guess starting parameters for, and its columns and range.

    Unlike ``SuggestPeakParamsRequest`` it names columns rather than carrying data: the Group Fit
    window holds no scan data, so ``GroupFitModel`` resolves it, the same way it does to fit.
    """

    source_scan_uuid: UUID
    y_col: str
    x_col: Optional[str] = None
    """``None`` uses the scan's own ``tavimeta.default_axis[0]``, as ``GroupFitRequest.x_col``."""
    range_min: str
    range_max: str
    """Raw panel text; a blank bound leaves that side of the scan open, as in a group fit."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


class SuggestGroupPeakParamsRequest(SuggestGroupParamsRequest):
    """Guess one peak's start for a scan of the group, for the panel's chosen peak ``shape``."""

    shape: str


class SuggestGroupBackgroundParamsRequest(SuggestGroupParamsRequest):
    """Guess the background's start for a scan of the group, for the panel's chosen ``background``."""

    background: str


class ColumnAverageSummary(BaseModel):
    """A column's per-scan mean and standard deviation (its spread during the scan)."""

    mean: float
    std: float


class GroupFitOutcome(BaseModel):
    """
    One scan of a group fit as just computed: its spec, and its result, curve, data and averages if it ran.

    ``result``/``curve``/``data`` are ``None`` when that scan was left out of the fit - no longer
    loaded, a bad range, a missing column, or a model that differs from the first scan's - so the
    user can see which scan to fix, as with ``FitOutcome``.
    """

    source_scan_uuid: UUID
    spec: FitSpec
    """The spec as requested, so the panel can show it back."""
    result: Optional[FitResultSummary] = None
    curve: Optional[FitCurve] = None
    data: Optional[FitData] = None
    """The (x, y, err) points this scan was fit against - its columns trimmed to its own fitting range."""
    averages: dict[str, ColumnAverageSummary] = {}
    """Averaged over the whole scan, not just its fitting range."""

    model_config = ConfigDict(arbitrary_types_allowed=True)
