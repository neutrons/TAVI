"""Fit request/result data model."""

from collections.abc import Container
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field

from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, UUIDFactory


class ParamField(BaseModel):
    """One editable parameter row's raw, unparsed text state: value / fixed / min / max / constraint."""

    value: str
    fixed: bool
    minimum: str
    maximum: str
    constrained: bool = False


class PeakField(BaseModel):
    """The peak box's raw field state, keyed to what ``Fit.fit`` needs (amplitude, center, width)."""

    shape: str
    amplitude: ParamField
    center: ParamField
    fwhm: ParamField


class FitSpec(BaseModel):
    """
    The fitting panel's raw field state: range, background, and peaks - how to fit, not what to fit.

    Carries no series, so one spec can be applied to every series in a batch; ``FitMember`` binds it
    to the one series it was actually fit against.
    """

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


class FitMember(FitSpec):
    """
    One series' share of a fit: the spec it was fit with, and the parameters that fit produced.

    ``series`` (rather than a bare source scan uuid + column names) is what lets a fit be
    recomputed later against whatever the source scan's data currently is - same pointer-not-data
    philosophy as ``Plot``/``PlotSeries`` itself, which never caches resolved arrays either.
    """

    series: PlotSeries
    result: Optional[FitResultSummary] = None
    """The fitted parameters from the last time this member was fit, or ``None`` if that fit failed.

    Only scalar parameters are kept - the use case asks for fitted parameters to be saved in the
    project - never the evaluated curve. This is what lets the panel show a member's results again
    when its series is picked from the "Current Plot" dropdown, without refitting."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @property
    def source_scan_uuid(self) -> UUID:
        """The scan this member was fit against - its identity within its ``FitEntry``."""
        return self.series.source_scan_uuid


class FitEntry(BaseModel):
    """
    A saved fit: one member per series it covers, in the order they were fit.

    A fit of one scan is simply a fit with one member; a sequential fit is several, each seeded from
    the one before. Once fit, nothing distinguishes the two - every member is recomputed from its own
    stored spec - so they share this one type, the same way a ``Plot`` of one series and a ``Plot``
    of many do.

    TAVI must not cache the actual fitted data points: selecting a fit later recomputes its
    curves fresh (see ``RecomputeFitEvent`` / ``FitModel``) against the source series' current data,
    rather than replaying stale points cached at fit time.
    """

    uuid: UUID = UUIDFactory()
    name: str = ""
    members: list[FitMember]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def member_for(self, source_scan_uuid: UUID) -> Optional[FitMember]:
        """Return the member fit against ``source_scan_uuid``, or ``None`` if this fit doesn't cover that scan."""
        return next((member for member in self.members if member.source_scan_uuid == source_scan_uuid), None)

    def with_members(self, members: list[FitMember]) -> "FitEntry":
        """Return a copy with each of ``members`` replacing this fit's member for the same scan, or appended."""
        merged = list(self.members)
        for member in members:
            index = next(
                (i for i, existing in enumerate(merged) if existing.source_scan_uuid == member.source_scan_uuid), None
            )
            if index is None:
                merged.append(member)
            else:
                merged[index] = member
        return self.model_copy(update={"members": merged})

    def without_scans(self, scan_uuids: Container[UUID]) -> Optional["FitEntry"]:
        """
        Return this fit with every member fit against ``scan_uuids`` dropped - mirrors ``Plot.without_scans``.

        Returns ``self`` when nothing matched, and ``None`` once no member is left.
        """
        surviving = [member for member in self.members if member.source_scan_uuid not in scan_uuids]
        if len(surviving) == len(self.members):
            return self
        if not surviving:
            return None
        return self.model_copy(update={"members": surviving})


class FitRequest(BaseModel):
    """
    Everything ``FitModel`` needs to fit one or more series: the panel's raw spec and the series to fit.

    Mirrors ``PlotFields``: the view hands over raw text, the model does all numeric
    parsing/validation and reports errors instead of failing silently. Data is not resolved here -
    ``FitModel`` resolves each series against its own scan handle, the same way it does to recompute.
    """

    spec: FitSpec
    series: list[PlotSeries]
    """Fit in this order. ``spec`` applies as-is to the first; see ``seed_from_previous`` for the rest."""
    seed_from_previous: bool = False
    """Start each series from the previous series' fitted values instead of ``spec``'s - a sequential
    fit. Only the values are carried over; fixed flags and bounds stay as ``spec`` has them."""
    fit_uuid: Optional[UUID] = None
    """The fit this request (re)fits, when the panel already has one for these series.

    ``None`` mints a fresh ``FitEntry``; supplying one instead refits in place - each fitted series
    replaces that fit's member for the same scan - so ``TaviProjectModel`` updates that entry rather
    than growing the project tree by one fit per click of Perform Fit."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


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


class FitCurve(BaseModel):
    """
    A freshly (re)computed member's evaluated curve, for plotting only - never cached.

    Published on ``SyncFitEvent``, whether the fit was just performed or recomputed after being
    selected from the project tree.
    """

    source_scan_uuid: UUID
    """The scan the fitted series came from. One drawn curve per source scan - the same key
    ``PlotterPresenter`` caches these under - so redrawing replaces that series' curve rather
    than stacking another one on top of it."""
    scan_name: str
    """Legend label for this curve - the fitted series' ``display_label``, so it names the run, not a shared title."""
    x: list[float]
    best_fit: list[float]
    components: dict[str, list[float]] = {}
    """Each model component (``bg_``, ``peak1_``, ...) evaluated on its own over the same ``x``,
    for "Plot Separately". Always carried, whether or not the panel is currently showing them, so
    toggling the checkbox redraws from what's already here instead of re-running the fit. A fit
    with a single component leaves this empty - that component *is* ``best_fit``."""

    model_config = ConfigDict(arbitrary_types_allowed=True)


class FitData(BaseModel):
    """The (x, y, err) points a member was fit against, resolved fresh from its series - for display only."""

    x: list[float]
    y: list[float]
    err: list[float]


class FitOutcome(BaseModel):
    """
    One member as just (re)computed: its spec and result, plus its curve and data if the fit ran.

    ``curve``/``data`` are ``None`` (and ``member.result`` too) when that member's fit raised -
    bad range, unsupported shape, a scan no longer loaded. The member is still kept, so the user
    can fix that one scan's parameters and refit it alone: the use case puts judging whether a fit
    worked on the user, not on TAVI.
    """

    member: FitMember
    curve: Optional[FitCurve] = None
    data: Optional[FitData] = None
