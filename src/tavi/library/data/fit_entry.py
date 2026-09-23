"""Fit request/result data model."""

from typing import Optional

from pydantic import BaseModel, ConfigDict, Field

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
    """The background line's constant term (its intercept)."""
    background_slope: ParamField = Field(
        default_factory=lambda: ParamField(value="", fixed=False, minimum="", maximum="")
    )
    """The background line's slope. Blank by default, which means "guess it from the data",
    exactly as a blank peak parameter does; fixing it at 0 gives a flat background."""
    peaks: list[PeakField]
    """One entry per peak panel, in panel order. Every entry becomes its own lmfit component
    (``peak1_``, ``peak2_``, ...) summed into one composite model."""

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
    fit_uuid: Optional[UUID] = None
    """The uuid of the fit this request re-runs, when the panel already has one for this series.

    ``None`` means "no fit yet" and mints a fresh ``FitEntry``; supplying one instead refits in
    place, so ``TaviProjectModel`` overwrites that entry rather than growing the project tree by
    one fit per click of Perform Fit."""


class SuggestPeakParamsRequest(BaseModel):
    """Everything ``FitModel`` needs to guess one peak's starting amplitude/center/FWHM from data."""

    source_scan_uuid: UUID
    x: list[float]
    y: list[float]
    range_min: str
    range_max: str
    shape: str

    model_config = ConfigDict(arbitrary_types_allowed=True)


class SuggestBackgroundParamsRequest(BaseModel):
    """
    Everything ``FitModel`` needs to guess the background's starting parameters from data.

    Mirrors ``SuggestPeakParamsRequest``; ``background`` is the raw combo text, resolved (and
    rejected if unsupported) by the model rather than the view.
    """

    source_scan_uuid: UUID
    x: list[float]
    y: list[float]
    range_min: str
    range_max: str
    background: str

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

    source_scan_uuid: UUID
    """The scan the fitted series came from. One drawn curve per source scan - the same key
    ``PlotterPresenter`` caches these under - so redrawing replaces that series' curve rather
    than stacking another one on top of it."""
    scan_name: str
    x: list[float]
    best_fit: list[float]
    components: dict[str, list[float]] = {}
    """Each model component (``bg_``, ``peak1_``, ...) evaluated on its own over the same ``x``,
    for "Plot Separately". Always carried, whether or not the panel is currently showing them, so
    toggling the checkbox redraws from what's already here instead of re-running the fit. A fit
    with a single component leaves this empty - that component *is* ``best_fit``."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


class PeakResult(BaseModel):
    """One fitted peak's values and 1-sigma uncertainties."""

    amplitude: float
    amplitude_err: Optional[float]
    center: float
    center_err: Optional[float]
    fwhm: float
    fwhm_err: Optional[float]


class FitResultSummary(BaseModel):
    """Scalar fit-result readback for the fitting panel: value + 1-sigma uncertainty per parameter."""

    reduced_chi_squared: float
    peaks: list[PeakResult]
    """One entry per requested peak, in the same order as ``FitSpec.peaks`` - so entry *i* reads
    back into peak panel *i*."""
    background_constant: Optional[float] = None
    background_constant_err: Optional[float] = None
    background_slope: Optional[float] = None
    background_slope_err: Optional[float] = None
