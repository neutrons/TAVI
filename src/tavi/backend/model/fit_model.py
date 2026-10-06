"""Model for computing 1D peak fits."""

import math
import threading
from collections.abc import Generator
from contextlib import contextmanager
from types import SimpleNamespace
from typing import Any, Optional

import numpy as np

from tavi.backend.model.interface.fit_model_interface import FitModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.library.data.fit_entry import (
    FitCurve,
    FitData,
    FitEntry,
    FitMember,
    FitOutcome,
    FitRequest,
    FitResultSummary,
    FitSpec,
    ParamField,
    PeakField,
    PeakResult,
    SuggestBackgroundParamsRequest,
    SuggestPeakParamsRequest,
)
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import UUID, RawScan, new_uuid
from tavi.library.fit import Fit, FitPackage, ModelName
from tavi.library.fit.fit import FitResult
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.event.type.presenter_event import (
    BackgroundParamsSuggestedEvent,
    FitRecomputeEvent,
    PeakParamsSuggestedEvent,
    SaveFitEvent,
    SyncFitEvent,
    SyncFitSpecEvent,
)
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError

# lmfit's Fit.fit expects a width in `sigma`, but the peak table shows FWHM. Gaussian's
# conversion is exact; Lorentzian's is exact for its own sigma definition; Voigt's true FWHM
# has no closed form (it depends on both its Gaussian and Lorentzian widths), so it reuses
# the Lorentzian approximation and leaves `gamma` at lmfit's default (tied to `sigma`).
_FWHM_TO_SIGMA = {
    ModelName.Gaussian: 1 / (2 * math.sqrt(2 * math.log(2))),
    ModelName.Lorentzian: 0.5,
    ModelName.Voigt: 0.5,
}

_PEAK_SHAPES = {
    "Gaussian": ModelName.Gaussian,
    "Lorentzian": ModelName.Lorentzian,
    "Voigt": ModelName.Voigt,
}

# The only background the panel offers. A flat background is this one with its slope fixed to
# zero, which the table's "fix" checkbox already expresses - no separate model needed.
_BACKGROUND = "Linear"


def _peak_prefix(index: int) -> str:
    """Return the lmfit prefix for peak ``index`` (0-based), numbered from 1 to match its panel label."""
    return f"peak{index + 1}_"


# A seed this far below the data's own scale is a fit that landed on (numerically) zero. lmfit's
# finite-difference step is proportional to the starting value, so seeding the next member with,
# say, 9e-14 leaves that parameter effectively frozen for the rest of the chain; starting from
# exactly 0 instead makes the fitter fall back to an absolute step.
_SEED_ZERO_FRACTION = 1e-6


def _seeded_value(field: ParamField, value: Optional[float], scale: float) -> ParamField:
    """Return ``field`` starting from ``value`` instead, keeping its fixed flag and bounds."""
    if value is None:
        return field
    if abs(value) < _SEED_ZERO_FRACTION * scale:
        value = 0.0
    return field.model_copy(update={"value": repr(float(value))})


def _seeded_spec(spec: FitSpec, result: FitResultSummary, data: FitData) -> FitSpec:
    """Return ``spec`` with every value replaced by ``result``'s - how one sequential member seeds the next."""
    x_scale = float(np.max(np.abs(data.x))) if data.x else 0.0
    x_span = float(np.ptp(data.x)) if data.x else 0.0
    y_scale = float(np.max(np.abs(data.y))) if data.y else 0.0
    slope_scale = y_scale / x_span if x_span else y_scale
    peaks = [
        peak.model_copy(
            update={
                "amplitude": _seeded_value(peak.amplitude, fitted.amplitude, y_scale * x_span),
                "center": _seeded_value(peak.center, fitted.center, x_scale),
                "fwhm": _seeded_value(peak.fwhm, fitted.fwhm, x_span),
            }
        )
        for peak, fitted in zip(spec.peaks, result.peaks)
    ]
    return spec.model_copy(
        update={
            "peaks": peaks,
            "background_constant": _seeded_value(spec.background_constant, result.background_constant, y_scale),
            "background_slope": _seeded_value(spec.background_slope, result.background_slope, slope_scale),
        }
    )


class FitModel(FitModelInterface):
    """Fits one or more series against a peak/background spec and announces the result."""

    def __init__(self, raw_scans: dict[UUID, RawScan], fits: dict[UUID, FitEntry]) -> None:
        """Init with live handles into TaviData's raw_scans (to resolve each series' data) and fits (to sync specs)."""
        self._raw_scans = raw_scans
        self._fits = fits
        # Set while fitting a batch, so each member's validation failures are gathered into one
        # report instead of popping one error dialog per scan. Thread-local because perform_fit
        # runs on the proxy's worker while a recompute runs on whichever thread published it.
        self._error_sink = threading.local()
        self._event_broker = EventBroker()
        self._event_broker.register(FitRecomputeEvent, self._handle_fit_recompute_event)

    def perform_fit(self, request: FitRequest) -> ModelResponse:
        """
        Fit each of ``request.series`` in order, then publish the members as one ``SaveFitEvent`` and ``SyncFitEvent``.

        With ``seed_from_previous`` set, every series after the first starts from the last
        successful fit's values - the sequential fit. A member whose fit fails is still published
        (with no result) so the user can fix and refit it alone, unless every member failed, in
        which case there's nothing to save and nothing is published.
        """
        spec: FitSpec = request.spec
        outcomes: list[FitOutcome] = []
        with self._collecting_errors(len(request.series) > 1) as errors:
            for series in request.series:
                member = FitMember(series=series, **spec.model_dump())
                errors.label = series.run_name
                outcome = self._compute_member(member)
                outcomes.append(outcome)
                if request.seed_from_previous and outcome.member.result is not None and outcome.data is not None:
                    spec = _seeded_spec(spec, outcome.member.result, outcome.data)

        if any(outcome.member.result is not None for outcome in outcomes):
            fit_uuid = request.fit_uuid if request.fit_uuid is not None else new_uuid()
            # Saved first, so the fit windows know the sync that follows is a fresh fit.
            self._event_broker.publish(
                SaveFitEvent(fit_uuid=fit_uuid, members=[outcome.member for outcome in outcomes])
            )
            self._event_broker.publish(SyncFitEvent(fit_uuid=fit_uuid, outcomes=outcomes))
        return ModelResponse(code=ResponseCode.OK)

    def sync_fit_spec(self, fit_uuid: UUID, source_scan_uuid: UUID) -> ModelResponse:
        """Publish one saved member's spec and result for the fitting panel - a miss means the fit or member is gone."""
        fit = self._fits.get(fit_uuid)
        member = fit.member_for(source_scan_uuid) if fit is not None else None
        if member is not None:
            self._event_broker.publish(SyncFitSpecEvent(fit_uuid=fit_uuid, member=member))
        return ModelResponse(code=ResponseCode.OK)

    def _handle_fit_recompute_event(self, e: FitRecomputeEvent) -> None:
        """Recompute every member of each selected fit from its own stored spec - never a cached curve, never re-seeded."""
        for fit in e.fits:
            members = [member for member in fit.members if member.source_scan_uuid in self._raw_scans]
            with self._collecting_errors(len(members) > 1) as errors:
                outcomes = []
                for member in members:
                    errors.label = member.series.run_name
                    outcomes.append(self._compute_member(member))
            if any(outcome.member.result is not None for outcome in outcomes):
                self._event_broker.publish(SyncFitEvent(fit_uuid=fit.uuid, outcomes=outcomes))

    def _compute_member(self, member: FitMember) -> FitOutcome:
        """Resolve one member's series, fit it, and return its outcome - with no result if the fit failed."""
        if member.source_scan_uuid not in self._raw_scans:
            self._report_error(f"Scan '{member.series.run_name}' is no longer loaded.")
            return FitOutcome(member=member.model_copy(update={"result": None}))
        x, y, err = resolve_series(member.series, self._raw_scans)
        computed = self._run_fit(member, x, y, err)
        if computed is None:
            return FitOutcome(member=member.model_copy(update={"result": None}))
        curve_x, result = computed

        # best_fit is evaluated only at the scan's own points, which draws as straight segments
        # on a coarse scan - re-evaluate on a fine grid, the same way browser.py plots a fit.
        x_fine = np.linspace(curve_x.min(), curve_x.max(), 300)
        curve = FitCurve(
            source_scan_uuid=member.source_scan_uuid,
            scan_name=member.series.display_label,
            x=x_fine.tolist(),
            best_fit=result.raw.eval(x=x_fine).tolist(),
            components=self._evaluate_components(result, x_fine),
        )
        summary = self._build_result_summary(member, result)
        return FitOutcome(
            member=member.model_copy(update={"result": summary}),
            curve=curve,
            data=FitData(x=x.tolist(), y=y.tolist(), err=err.tolist()),
        )

    @contextmanager
    def _collecting_errors(self, enabled: bool) -> Generator[Any, None, None]:
        """Gather ``_report_error`` messages, each prefixed by ``sink.label``, into one report when ``enabled``."""
        if not enabled:
            yield SimpleNamespace(label="")
            return
        sink = SimpleNamespace(label="", messages=[])
        self._error_sink.current = sink
        try:
            yield sink
        finally:
            self._error_sink.current = None
            if sink.messages:
                self._publish_error("\n".join(sink.messages))

    def _evaluate_components(self, result: FitResult, x_fine: np.ndarray) -> dict[str, list[float]]:
        """Evaluate each model component separately over ``x_fine``, the way browser.py's show_components does."""
        # A single-component fit has nothing to separate out - its one component is the composite
        # curve already being published, and drawing it twice would just overplot it.
        if len(result.components) < 2:
            return {}
        return {prefix: np.asarray(values).tolist() for prefix, values in result.raw.eval_components(x=x_fine).items()}

    def _build_result_summary(self, spec: FitSpec, result: FitResult) -> FitResultSummary:
        """Read every fitted peak (and, if present, the background) out of a FitResult."""
        peaks = []
        for index in range(len(spec.peaks)):
            peak = result[_peak_prefix(index)]
            peaks.append(
                PeakResult(
                    amplitude=peak["amplitude"],
                    amplitude_err=peak.errors["amplitude"],
                    center=peak["center"],
                    center_err=peak.errors["center"],
                    fwhm=peak["fwhm"],
                    fwhm_err=peak.errors["fwhm"],
                )
            )
        summary_kwargs: dict[str, Any] = dict(reduced_chi_squared=result.reduced_chi_squared, peaks=peaks)
        if spec.background == _BACKGROUND:
            background = result["bg_"]
            summary_kwargs["background_constant"] = background["intercept"]
            summary_kwargs["background_constant_err"] = background.errors["intercept"]
            summary_kwargs["background_slope"] = background["slope"]
            summary_kwargs["background_slope_err"] = background.errors["slope"]
        return FitResultSummary(**summary_kwargs)

    def suggest_peak_params(self, request: SuggestPeakParamsRequest) -> ModelResponse:
        """Trim to the fitting range, guess amplitude/center/FWHM from data, and publish the result."""
        x = np.array(request.x)
        y = np.array(request.y)
        mask = self._mask_range(request.range_min, request.range_max, x)
        if mask is None:
            return ModelResponse(code=ResponseCode.OK)
        x, y = x[mask], y[mask]

        shape = self._resolve_peak_shape(request.shape)
        if shape is None:
            return ModelResponse(code=ResponseCode.OK)

        guess = Fit(FitPackage.lmfit).guess(x, y, shape, prefix="peak_")
        self._event_broker.publish(
            PeakParamsSuggestedEvent(
                source_scan_uuid=request.source_scan_uuid,
                amplitude=guess["amplitude"],
                center=guess["center"],
                fwhm=guess["fwhm"],
            )
        )
        return ModelResponse(code=ResponseCode.OK)

    def suggest_background_params(self, request: SuggestBackgroundParamsRequest) -> ModelResponse:
        """Trim to the fitting range, guess the background's slope/intercept, and publish the result."""
        x = np.array(request.x)
        y = np.array(request.y)
        mask = self._mask_range(request.range_min, request.range_max, x)
        if mask is None:
            return ModelResponse(code=ResponseCode.OK)
        x, y = x[mask], y[mask]

        if request.background != _BACKGROUND:
            self._report_error(f"Background '{request.background}' is not supported yet.")
            return ModelResponse(code=ResponseCode.OK)

        # lmfit's linear guess is a least-squares line through every point in range, peak
        # included, so it sits high when the peak is strong. That is fine for a starting
        # value - the fit refines it - but it is not a peak-free background estimate.
        guess = Fit(FitPackage.lmfit).guess(x, y, ModelName.Linear, prefix="bg_")
        self._event_broker.publish(
            BackgroundParamsSuggestedEvent(
                source_scan_uuid=request.source_scan_uuid,
                slope=guess["slope"],
                intercept=guess["intercept"],
            )
        )
        return ModelResponse(code=ResponseCode.OK)

    def _run_fit(
        self, spec: FitSpec, x: np.ndarray, y: np.ndarray, err: np.ndarray
    ) -> Optional[tuple[np.ndarray, FitResult]]:
        """Trim to the fitting range, build the composite model, and run a weighted fit."""
        mask = self._mask_range(spec.range_min, spec.range_max, x)
        if mask is None:
            return None
        x, y, err = x[mask], y[mask], err[mask]

        model_dict = self._build_model_dict(spec, x, y)
        if model_dict is None:
            return None

        try:
            result = Fit(FitPackage.lmfit).fit(x=x, y=y, model_dict=model_dict, err=err)
        except ValueError as error:
            self._report_error(str(error))
            return None
        return x, result

    def _mask_range(self, range_min_text: str, range_max_text: str, x: np.ndarray) -> Optional[np.ndarray]:
        """Parse the fitting range and return a boolean mask into ``x``, or report an error and return None."""
        try:
            range_min = float(range_min_text)
            range_max = float(range_max_text)
        except ValueError:
            self._report_error(f"Fitting range '{range_min_text}' to '{range_max_text}' is not numeric.")
            return None

        mask = (x >= range_min) & (x <= range_max)
        if mask.sum() < 2:
            self._report_error(f"Fitting range {range_min} to {range_max} contains too few points to fit.")
            return None
        return mask

    def _build_model_dict(
        self, spec: FitSpec, x: np.ndarray, y: np.ndarray
    ) -> Optional[list[tuple[ModelName, dict[str, Any]]]]:
        """Translate the spec's background/peak fields into a Fit.fit model_dict, or None on error."""
        model_dict: list[tuple[ModelName, dict[str, Any]]] = []

        if spec.background not in ("None", _BACKGROUND):
            self._report_error(f"Background '{spec.background}' is not supported yet.")
            return None
        if spec.background == _BACKGROUND:
            background_component = self._build_background_linear(spec, x, y)
            if background_component is None:
                return None
            model_dict.append(background_component)

        if not spec.peaks:
            self._report_error("A fit needs at least one peak.")
            return None
        for index, peak in enumerate(spec.peaks):
            peak_component = self._build_peak(peak, _peak_prefix(index), x, y)
            if peak_component is None:
                return None
            model_dict.append(peak_component)

        return model_dict

    def _build_background_linear(
        self, spec: FitSpec, x: np.ndarray, y: np.ndarray
    ) -> Optional[tuple[ModelName, dict[str, Any]]]:
        """Return the "Linear" background with both terms free, or None on a parse error."""
        slope = self._parse_param("background slope", spec.background_slope)
        intercept = self._parse_param("background constant", spec.background_constant)
        if slope is None or intercept is None:
            return None
        slope_value, slope_min, slope_max = slope
        intercept_value, intercept_min, intercept_max = intercept

        if slope_value is None or intercept_value is None:
            # At least one was left blank - guess every unset one from the data itself, the
            # same fallback _build_peak applies to a blank amplitude/center/FWHM.
            guess = Fit(FitPackage.lmfit).guess(x, y, ModelName.Linear, prefix="bg_")
            slope_value = guess["slope"] if slope_value is None else slope_value
            intercept_value = guess["intercept"] if intercept_value is None else intercept_value

        return (
            ModelName.Linear,
            {
                "prefix": "bg_",
                "slope": slope_value,
                "intercept": intercept_value,
                "set": {
                    "slope": self._bounds(not spec.background_slope.fixed, slope_min, slope_max),
                    "intercept": self._bounds(not spec.background_constant.fixed, intercept_min, intercept_max),
                },
            },
        )

    def _build_peak(
        self, peak: PeakField, prefix: str, x: np.ndarray, y: np.ndarray
    ) -> Optional[tuple[ModelName, dict[str, Any]]]:
        """
        Return one peak's (ModelName, params) component under ``prefix``, or None on error.

        Note that a blank field falls back to a guess made against the *whole* range, which is
        the same guess for every peak - so leaving several peaks entirely blank starts them all
        on top of each other. Give each one a center (or use its Suggest Params. button).
        """
        shape = self._resolve_peak_shape(peak.shape)
        if shape is None:
            return None

        amplitude = self._parse_param("amplitude", peak.amplitude)
        center = self._parse_param("center", peak.center)
        fwhm = self._parse_param("FWHM", peak.fwhm)
        if amplitude is None or center is None or fwhm is None:
            return None

        to_sigma = _FWHM_TO_SIGMA[shape]
        amplitude_value, amplitude_min, amplitude_max = amplitude
        center_value, center_min, center_max = center
        fwhm_value, fwhm_min, fwhm_max = fwhm

        if amplitude_value is None or center_value is None or fwhm_value is None:
            # At least one was left blank - guess every unset one from the data itself.
            guess = Fit(FitPackage.lmfit).guess(x, y, shape, prefix=prefix)
            amplitude_value = guess["amplitude"] if amplitude_value is None else amplitude_value
            center_value = guess["center"] if center_value is None else center_value
            fwhm_value = guess["fwhm"] if fwhm_value is None else fwhm_value

        sigma_value = fwhm_value * to_sigma
        sigma_min = fwhm_min * to_sigma if fwhm_min is not None else None
        sigma_max = fwhm_max * to_sigma if fwhm_max is not None else None

        return (
            shape,
            {
                "prefix": prefix,
                "amplitude": amplitude_value,
                "center": center_value,
                "sigma": sigma_value,
                "set": {
                    "amplitude": self._bounds(not peak.amplitude.fixed, amplitude_min, amplitude_max),
                    "center": self._bounds(not peak.center.fixed, center_min, center_max),
                    "sigma": self._bounds(not peak.fwhm.fixed, sigma_min, sigma_max),
                },
            },
        )

    def _resolve_peak_shape(self, shape_text: str) -> Optional[ModelName]:
        """Map a peak_shape_combo text to its ModelName, reporting an error for an unsupported one."""
        shape = _PEAK_SHAPES.get(shape_text)
        if shape is None:
            self._report_error(f"Peak shape '{shape_text}' is not supported yet.")
        return shape

    def _bounds(self, vary: bool, minimum: Optional[float], maximum: Optional[float]) -> dict[str, Any]:
        """Build an lmfit ``Parameter.set`` options dict from a vary flag and optional bounds."""
        options: dict[str, Any] = {"vary": vary}
        if minimum is not None:
            options["min"] = minimum
        if maximum is not None:
            options["max"] = maximum
        return options

    def _parse_param(
        self, label: str, field: ParamField
    ) -> Optional[tuple[Optional[float], Optional[float], Optional[float]]]:
        """
        Parse a ParamField's value/min/max text into floats, reporting an error and returning None on a bad parse.

        A blank ``value`` isn't an error: it returns ``None`` for the value (min/max still parsed)
        so the caller can fall back to a heuristic guess - "TAS user does not set initial fitting
        parameters" is meant to be supported, not rejected.
        """
        try:
            value = float(field.value) if field.value.strip() else None
            minimum = float(field.minimum) if field.minimum.strip() else None
            maximum = float(field.maximum) if field.maximum.strip() else None
        except ValueError:
            self._report_error(f"{label} value/min/max must be numeric.")
            return None
        return value, minimum, maximum

    def _report_error(self, message: str) -> None:
        """Surface a fit validation failure to the user instead of failing silently."""
        sink = getattr(self._error_sink, "current", None)
        if sink is not None:
            sink.messages.append(f"{sink.label}: {message}")
            return
        self._publish_error(message)

    def _publish_error(self, message: str) -> None:
        self._event_broker.publish(ExceptionEvent(error=NonRecoverableError(message, "")))
