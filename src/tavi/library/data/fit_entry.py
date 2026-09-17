"""Fit request/result data model."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, UUIDFactory


class ParamField(BaseModel):
    """One editable parameter row's raw, unparsed text state: value / fixed / min / max."""

    value: str
    fixed: bool
    minimum: str
    maximum: str


class PeakField(BaseModel):
    """The peak box's raw field state, keyed to what ``Fit.fit`` needs (amplitude, center, width)."""

    shape: str
    amplitude: ParamField
    center: ParamField
    fwhm: ParamField


class FitSpec(BaseModel):
    """
    Everything needed to redo a fit against a series' live data, short of the resolved (x, y, err) arrays.

    ``series`` (rather than a bare source scan uuid + column names) is what lets a fit be
    recomputed later against whatever the source scan's data currently is - same pointer-not-data
    philosophy as ``Plot``/``PlotSeries`` itself, which never caches resolved arrays either.
    """

    series: PlotSeries
    range_min: str
    range_max: str
    background: str
    background_constant: ParamField
    peak: PeakField

    model_config = ConfigDict(arbitrary_types_allowed=True)


class FitRequest(FitSpec):
    """
    Everything ``FitModel`` needs to run one fit: already-resolved data plus raw (unparsed) field state.

    Mirrors ``PlotFields``: the view hands over raw text, the model does all numeric
    parsing/validation and reports errors instead of failing silently.
    """

    x: list[float]
    y: list[float]
    err: list[float]


class SuggestPeakParamsRequest(BaseModel):
    """Everything ``FitModel`` needs to guess one peak's starting amplitude/center/FWHM from data."""

    source_scan_uuid: UUID
    x: list[float]
    y: list[float]
    range_min: str
    range_max: str
    shape: str

    model_config = ConfigDict(arbitrary_types_allowed=True)


class FitEntry(FitSpec):
    """
    A saved fit's specification - what to fit and how, not its evaluated curve.

    TAVI must not cache the actual fitted data points: selecting a fit later recomputes its
    curve fresh (see ``FitRecomputeEvent`` / ``FitModel``) against the source series' current data,
    rather than replaying stale points cached at fit time.
    """

    uuid: UUID = UUIDFactory()


class FitCurve(BaseModel):
    """
    A freshly (re)computed fit's evaluated curve, for plotting only - never cached.

    Published alongside a ``FitEntry`` on ``FitComputedEvent``, whether the fit was just performed
    or recomputed after being selected from the project tree.
    """

    scan_name: str
    x: list[float]
    best_fit: list[float]

    model_config = ConfigDict(arbitrary_types_allowed=True)


class FitResultSummary(BaseModel):
    """Scalar fit-result readback for the fitting panel: value + 1-sigma uncertainty per parameter."""

    reduced_chi_squared: float
    amplitude: float
    amplitude_err: Optional[float]
    center: float
    center_err: Optional[float]
    fwhm: float
    fwhm_err: Optional[float]
    background_constant: Optional[float] = None
    background_constant_err: Optional[float] = None
