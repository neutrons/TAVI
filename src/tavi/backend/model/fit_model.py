"""Model for computing 1D peak fits."""

import math
from typing import Any, Optional

import numpy as np

from tavi.backend.model.interface.fit_model_interface import FitModelInterface
from tavi.backend.model.plot_resolver import resolve_series
from tavi.library.data.fit_entry import (
    FitCurve,
    FitEntry,
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
from tavi.library.data.scan import UUID, RawScan
from tavi.library.fit import Fit, FitPackage, ModelName
from tavi.library.fit.fit import FitResult
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.event.type.presenter_event import (
    BackgroundParamsSuggestedEvent,
    FitComputedEvent,
    FitRecomputeEvent,
    PeakParamsSuggestedEvent,
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


class FitModel(FitModelInterface):
    """Computes a peak fit from an already-resolved (x, y, err) series and announces the result."""

    def __init__(self, raw_scans: dict[UUID, RawScan]) -> None:
        """Init with a live handle into TaviData's raw_scans storage, for recomputing a selected fit."""
        self._raw_scans = raw_scans
        self._event_broker = EventBroker()
        self._event_broker.register(FitRecomputeEvent, self._handle_fit_recompute_event)

    def perform_fit(self, request: FitRequest) -> ModelResponse:
        """Run a fresh fit from an already-resolved request and publish a FitComputedEvent on success."""
        x = np.array(request.x)
        y = np.array(request.y)
        err = np.array(request.err)

        computed = self._run_fit(request, x, y, err)
        if computed is None:
            return ModelResponse(code=ResponseCode.OK)
        curve_x, result = computed

        fit_entry = FitEntry(
            series=request.series,
            range_min=request.range_min,
            range_max=request.range_max,
            background=request.background,
            background_constant=request.background_constant,
            background_slope=request.background_slope,
            peaks=request.peaks,
        )
        if request.fit_uuid is not None:
            # Refitting a series the panel already has a fit for: keep that uuid so the new spec
            # replaces it in TaviData.fits, instead of each click leaving another fit behind.
            fit_entry = fit_entry.model_copy(update={"uuid": request.fit_uuid})
        self._publish_fit_computed(fit_entry, curve_x, result)
        return ModelResponse(code=ResponseCode.OK)

    def _handle_fit_recompute_event(self, e: FitRecomputeEvent) -> None:
        """Recompute each selected fit against its source series' current data - never against a cached curve."""
        for fit in e.fits:
            if fit.series.source_scan_uuid not in self._raw_scans:
                continue
            x, y, err = resolve_series(fit.series, self._raw_scans)
            computed = self._run_fit(fit, x, y, err)
            if computed is None:
                continue
            curve_x, result = computed
            self._publish_fit_computed(fit, curve_x, result)

    def _publish_fit_computed(self, fit_entry: FitEntry, curve_x: np.ndarray, result: FitResult) -> None:
        """Build the transient curve/result payloads and publish them alongside the (cacheable) spec."""
        # best_fit is evaluated only at the scan's own points, which draws as straight segments
        # on a coarse scan - re-evaluate on a fine grid, the same way browser.py plots a fit.
        x_fine = np.linspace(curve_x.min(), curve_x.max(), 300)
        curve = FitCurve(
            source_scan_uuid=fit_entry.series.source_scan_uuid,
            scan_name=fit_entry.series.scan_name,
            x=x_fine.tolist(),
            best_fit=result.raw.eval(x=x_fine).tolist(),
        )
        summary = self._build_result_summary(fit_entry, result)
        self._event_broker.publish(FitComputedEvent(fit=fit_entry, curve=curve, result=summary))

    def _build_result_summary(self, fit_entry: FitEntry, result: FitResult) -> FitResultSummary:
        """Read every fitted peak (and, if present, the background) out of a FitResult."""
        peaks = []
        for index in range(len(fit_entry.peaks)):
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
        if fit_entry.background == _BACKGROUND:
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
        self._event_broker.publish(ExceptionEvent(error=NonRecoverableError(message, "")))
