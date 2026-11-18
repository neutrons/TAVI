"""Model for fitting a group of scans: each scan's panel spec, run through one ``GroupFit``."""

import threading
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from tavi.backend.model.fit_model import _BACKGROUND, _COSINE_PARAMS, FitModel, _peak_prefix
from tavi.backend.model.interface.group_fit_model_interface import GroupFitModelInterface
from tavi.library.data.fit_entry import (
    FitCurve,
    FitData,
    FitSpec,
    FitSummary,
    FitSummaryPoint,
    ParamField,
    SuggestBackgroundParamsRequest,
    SuggestPeakParamsRequest,
)
from tavi.library.data.group_fit_entry import (
    ColumnAverageSummary,
    GroupFitOutcome,
    GroupFitRequest,
    SuggestGroupBackgroundParamsRequest,
    SuggestGroupParamsRequest,
    SuggestGroupPeakParamsRequest,
    SummarizeGroupFitRequest,
)
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import UUID, Scan
from tavi.library.fit.group_fit import GroupFit, GroupFitParams, GroupFitResult, ModelDict
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.event_interface import Event
from tavi.meta.event.type.presenter_event import (
    SaveFitSummaryEvent,
    SyncBackgroundParamsEvent,
    SyncCosineParamsEvent,
    SyncGroupBackgroundParamsEvent,
    SyncGroupCosineParamsEvent,
    SyncGroupFitEvent,
    SyncGroupPeakParamsEvent,
    SyncPeakParamsEvent,
)


class _GroupGuessBroker:
    """
    The translator's broker: forwards what it publishes, renaming single-fit guesses to their group-fit events.

    ``FitModel.suggest_*`` publish ``Sync*ParamsEvent``, which the 1D Fitting tab fills its panel
    from. A guess for the Group Fit window must reach that window only, so it goes out under the
    group's own event type instead. Everything else - error reports - passes through unchanged.
    """

    _GROUP_TYPES: dict[type, type] = {
        SyncPeakParamsEvent: SyncGroupPeakParamsEvent,
        SyncCosineParamsEvent: SyncGroupCosineParamsEvent,
        SyncBackgroundParamsEvent: SyncGroupBackgroundParamsEvent,
    }

    def publish(self, event: Event) -> None:
        """Publish ``event`` on the application's broker, as its group-fit type if it is a guess."""
        group_type = self._GROUP_TYPES.get(type(event))
        EventBroker().publish(group_type(**event.model_dump()) if group_type is not None else event)


class _FitSpecTranslator(FitModel):
    """
    ``FitModel``'s spec -> lmfit translation and result readback, without ``FitModel``'s event subscription.

    Borrowed by inheritance so a group fit parses a panel spec exactly as a single fit does - range
    mask, blank fields guessed from the data, FWHM -> sigma, fixed/bounds - without a second copy
    of that code to drift.
    """

    def __init__(self) -> None:
        """Set only what the borrowed helpers read."""
        # super().__init__() is skipped on purpose: a FitModel registers for RecomputeFitEvent, and a
        # second one would recompute and sync every saved fit a second time.
        self._raw_scans = {}
        self._fits = {}
        self._error_sink = threading.local()
        self._event_broker = _GroupGuessBroker()


@dataclass
class _PreparedScan:
    """One scan ready to fit: its spec, its x/y trimmed to that spec's range, and its model built against them."""

    uuid: UUID
    name: str
    spec: FitSpec
    x: np.ndarray
    y: np.ndarray
    model_dict: ModelDict
    trimmed_scan: Scan
    """The scan with only its x and y columns trimmed - what GroupFit reads, since it has no fitting range."""


def _open_range(text: str, unbounded: str) -> str:
    """Return a range bound's text, or ``unbounded`` when blank - a group's scans each span their own x."""
    return text.strip() or unbounded


def _column_key(scan: Scan, name: str) -> str:
    """Return the key ``name`` is stored under, by ``GroupFit._column``'s rule, so a trimmed copy replaces that column."""
    data = scan.data.data
    if name not in data and name[:1].isdigit() and f"_{name}" in data:
        return f"_{name}"
    return name


def _structure(model_dict: ModelDict) -> list[tuple[Any, str]]:
    """Return each component's (model, prefix) - what a group's scans must share for one model to cover them."""
    return [(model_name, params.get("prefix", "")) for model_name, params in model_dict]


def _param_slots(spec: FitSpec) -> list[tuple[ParamField, str, str]]:
    """Return every ParamField ``spec`` fits, with the (component prefix, model parameter) it becomes."""
    slots: list[tuple[ParamField, str, str]] = []
    if spec.background == _BACKGROUND:
        slots += [(spec.background_slope, "bg_", "slope"), (spec.background_constant, "bg_", "intercept")]
    for index, peak in enumerate(spec.peaks):
        prefix = _peak_prefix(index)
        if peak.cosine is not None:
            slots += [(getattr(peak.cosine, name), prefix, name) for name in _COSINE_PARAMS]
        elif peak.amplitude is not None and peak.center is not None and peak.fwhm is not None:
            # The table edits FWHM, but FitModel hands lmfit a sigma - the override must name what lmfit fits.
            slots += [
                (peak.amplitude, prefix, "amplitude"),
                (peak.center, prefix, "center"),
                (peak.fwhm, prefix, "sigma"),
            ]
    return slots


def _full_patches(model_dict: ModelDict) -> dict[str, dict[str, Any]]:
    """Return every component of ``model_dict`` as an override, so a scan is fit from its own spec alone."""
    return {
        params.get("prefix", ""): {key: value for key, value in params.items() if key != "prefix"}
        for _model_name, params in model_dict
    }


def _changed_patches(group_spec: FitSpec, scan: _PreparedScan) -> dict[str, dict[str, Any]]:
    """
    Return overrides for only the parameters ``scan``'s spec changes from ``group_spec``.

    A changed parameter overrides both its start and its fixed/bounds, read from the scan's own
    model so FWHM -> sigma and the constraint checkbox are applied as a single fit would. A changed
    parameter whose value was left blank overrides only fixed/bounds, leaving the chained seed to
    supply its start.
    """
    components = {params.get("prefix", ""): params for _model_name, params in scan.model_dict}
    patches: dict[str, dict[str, Any]] = {}
    # zip lines the slots up one to one: scans reach here only with the first scan's structure.
    for (group_field, _, _), (field, prefix, key) in zip(_param_slots(group_spec), _param_slots(scan.spec)):
        if field == group_field:
            continue
        component = components[prefix]
        patch = patches.setdefault(prefix, {})
        patch.setdefault("set", {})[key] = component["set"][key]
        if field.value.strip():
            patch[key] = component[key]
    return patches


def _group_fit_params(request: GroupFitRequest, prepared: list[_PreparedScan]) -> GroupFitParams:
    """
    Translate a request into ``GroupFitParams``: the first scan's model for the group, the rest as overrides.

    With ``chain`` on, a later scan overrides only what the user changed for it, so its other
    parameters start from the previous scan's fit. With it off, each scan overrides everything and
    is fit from its own spec.
    """
    first = prepared[0]
    overrides: dict[UUID, dict[str, dict[str, dict[str, Any]]]] = {}
    for scan in prepared[1:]:
        patches = _changed_patches(first.spec, scan) if request.chain else _full_patches(scan.model_dict)
        if patches:
            overrides[scan.uuid] = {request.y_col: patches}
    return GroupFitParams(
        x_col=request.x_col,
        fits={request.y_col: first.model_dict},
        averages=list(request.averages),
        overrides=overrides,
        chain=request.chain,
    )


class GroupFitModel(GroupFitModelInterface):
    """Fits a group of scans, each from its own panel spec, through ``GroupFit``, and announces the results."""

    def __init__(self, raw_scans: dict[UUID, Scan], processed_scans: dict[UUID, Scan]) -> None:
        """Init with live handles into TaviData's raw and processed scans, so combined scans can be grouped too."""
        self._raw_scans = raw_scans
        self._processed_scans = processed_scans
        self._translator = _FitSpecTranslator()
        self._event_broker = EventBroker()

    def perform_group_fit(self, request: GroupFitRequest) -> ModelResponse:
        """
        Fit every scan in ``request``, in order, and publish one ``SyncGroupFitEvent`` if any of them was fit.

        A scan that can't be fit - not loaded, bad range or column, or a model unlike the first
        scan's - is reported and left out, but still listed without a result. Every problem in one
        run is gathered into a single error report.
        """
        translator = self._translator
        if not request.y_col.strip():
            translator._report_error("Set the y column to fit.")
            return ModelResponse(code=ResponseCode.OK)
        with translator._collecting_errors(True) as errors:
            prepared: list[_PreparedScan] = []
            for uuid, spec in request.specs.items():
                errors.label = uuid.value
                scan = self._prepare(uuid, spec, request, errors)
                if scan is None:
                    continue
                if prepared and _structure(scan.model_dict) != _structure(prepared[0].model_dict):
                    translator._report_error(
                        f"Its background and peaks differ from {prepared[0].name}'s; a group fit shares one model."
                    )
                    continue
                prepared.append(scan)

            results: dict[UUID, GroupFitResult] = {}
            if prepared:
                errors.label = "Group fit"
                try:
                    fitted = GroupFit([scan.trimmed_scan for scan in prepared], _group_fit_params(request, prepared))
                    results = {result.uuid: result for result in fitted.fit()}
                except (KeyError, ValueError) as error:
                    # GroupFit fits the whole group in one call, so one failure leaves no result to show.
                    translator._report_error(str(error.args[0]) if error.args else str(error))
            elif not request.specs:
                translator._report_error("A group fit needs at least one scan.")

        if not results:
            return ModelResponse(code=ResponseCode.OK)

        by_uuid = {scan.uuid: scan for scan in prepared}
        outcomes = [
            self._outcome(by_uuid[uuid], results[uuid], request.y_col)
            if uuid in results
            else GroupFitOutcome(source_scan_uuid=uuid, spec=spec)
            for uuid, spec in request.specs.items()
        ]
        self._event_broker.publish(SyncGroupFitEvent(x_col=request.x_col, y_col=request.y_col, outcomes=outcomes))
        return ModelResponse(code=ResponseCode.OK)

    def summarize_group_fit(self, request: SummarizeGroupFitRequest) -> ModelResponse:
        """
        Pair each scan's fitted ``parameter`` with its ``x_col`` average and publish them as one ``SaveFitSummaryEvent``.

        A scan that can't be summarized - not loaded, no such column, or no such parameter - is
        reported and left out. Every problem in one run is gathered into a single error report.
        """
        translator = self._translator
        x_col = request.x_col.strip()
        if not x_col:
            translator._report_error("Set the column to average.")
            return ModelResponse(code=ResponseCode.OK)
        if not request.results:
            translator._report_error("Run a group fit before summarizing it.")
            return ModelResponse(code=ResponseCode.OK)

        points: list[FitSummaryPoint] = []
        with translator._collecting_errors(True) as errors:
            for uuid, result in request.results.items():
                errors.label = uuid.value
                scan = self._find_scan(uuid)
                if scan is None:
                    translator._report_error("Scan is no longer loaded.")
                    continue
                name = scan.tavimeta.friendly_name
                errors.label = name
                key = _column_key(scan, x_col)
                if key not in scan.data.data:
                    translator._report_error(f"Column {x_col!r} not found; available: {list(scan.data.data)}.")
                    continue
                parameter = result.parameters().get(request.parameter)
                if parameter is None:
                    translator._report_error(f"Its fit has no {request.parameter!r}.")
                    continue
                values = np.asarray(scan.data.data[key], dtype=float)
                y, y_err = parameter
                # The same mean and spread GroupFit gives an averaged column, over the whole scan.
                points.append(
                    FitSummaryPoint(
                        source_scan_uuid=uuid,
                        scan_name=name,
                        x=float(np.mean(values)),
                        x_err=float(np.std(values)),
                        y=y,
                        y_err=y_err,
                    )
                )

        if points:
            summary = FitSummary(
                name=f"{request.parameter} vs {x_col}", x_label=x_col, y_label=request.parameter, points=points
            )
            self._event_broker.publish(SaveFitSummaryEvent(summary=summary))
        return ModelResponse(code=ResponseCode.OK)

    def _find_scan(self, uuid: UUID) -> Optional[Scan]:
        if uuid in self._raw_scans:
            return self._raw_scans[uuid]
        return self._processed_scans.get(uuid)

    def suggest_peak_params(self, request: SuggestGroupPeakParamsRequest) -> ModelResponse:
        """Guess one peak's start from the scan's own columns, through ``FitModel.suggest_peak_params``."""
        resolved = self._resolve_for_guess(request)
        if resolved is not None:
            x, y = resolved
            self._translator.suggest_peak_params(
                SuggestPeakParamsRequest(shape=request.shape, **self._guess_fields(request, x, y))
            )
        return ModelResponse(code=ResponseCode.OK)

    def suggest_background_params(self, request: SuggestGroupBackgroundParamsRequest) -> ModelResponse:
        """Guess the background's start from the scan's own columns, through ``FitModel.suggest_background_params``."""
        resolved = self._resolve_for_guess(request)
        if resolved is not None:
            x, y = resolved
            self._translator.suggest_background_params(
                SuggestBackgroundParamsRequest(background=request.background, **self._guess_fields(request, x, y))
            )
        return ModelResponse(code=ResponseCode.OK)

    def _resolve_for_guess(self, request: SuggestGroupParamsRequest) -> Optional[tuple[np.ndarray, np.ndarray]]:
        """Resolve the scan a guess is for to its x/y, or report why it can't be."""
        if not request.y_col.strip():
            self._translator._report_error("Set the y column to fit.")
            return None
        scan = self._find_scan(request.source_scan_uuid)
        if scan is None:
            self._translator._report_error("Scan is no longer loaded.")
            return None
        resolved = self._columns(scan, request.x_col, request.y_col)
        return None if resolved is None else resolved[2:]

    @staticmethod
    def _guess_fields(request: SuggestGroupParamsRequest, x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
        """The fields a single-fit suggest request shares with a group one, its range opened where blank."""
        return dict(
            source_scan_uuid=request.source_scan_uuid,
            x=x.tolist(),
            y=y.tolist(),
            range_min=_open_range(request.range_min, "-inf"),
            range_max=_open_range(request.range_max, "inf"),
        )

    def _columns(
        self, scan: Scan, x_col: Optional[str], y_col: str
    ) -> Optional[tuple[str, str, np.ndarray, np.ndarray]]:
        """Return ``scan``'s (x key, y key, x, y) - ``x_col`` None meaning its default axis - or report a missing one."""
        data = scan.data.data
        x_key = _column_key(scan, x_col or scan.tavimeta.default_axis[0])
        y_key = _column_key(scan, y_col)
        missing = [key for key in (x_key, y_key) if key not in data]
        if missing:
            self._translator._report_error(f"Column {missing[0]!r} not found; available: {list(data)}.")
            return None
        return x_key, y_key, np.asarray(data[x_key], dtype=float), np.asarray(data[y_key], dtype=float)

    def _prepare(self, uuid: UUID, spec: FitSpec, request: GroupFitRequest, errors: Any) -> Optional[_PreparedScan]:
        """Resolve one scan's columns, trim them to its range, and build its model - or report why not."""
        translator = self._translator
        scan = self._find_scan(uuid)
        if scan is None:
            translator._report_error("Scan is no longer loaded.")
            return None
        name = scan.tavimeta.friendly_name
        errors.label = name

        data = scan.data.data
        resolved = self._columns(scan, request.x_col, request.y_col)
        if resolved is None:
            return None
        x_key, y_key, x, y = resolved
        # A blank bound is the whole scan: a group's scans each span their own x, so the Group Fit
        # window leaves the range open rather than guess one for all of them.
        mask = translator._mask_range(_open_range(spec.range_min, "-inf"), _open_range(spec.range_max, "inf"), x)
        if mask is None:
            return None
        x, y = x[mask], y[mask]
        model_dict = translator._build_model_dict(spec, x, y)
        if model_dict is None:
            return None

        # Only x and y are trimmed, so the averaged columns still cover the whole scan.
        trimmed_data = scan.data.model_copy(update={"data": {**data, x_key: x.tolist(), y_key: y.tolist()}})
        return _PreparedScan(
            uuid=uuid,
            name=name,
            spec=spec,
            x=x,
            y=y,
            model_dict=model_dict,
            trimmed_scan=scan.model_copy(update={"data": trimmed_data}),
        )

    def _outcome(self, scan: _PreparedScan, result: GroupFitResult, y_col: str) -> GroupFitOutcome:
        """Read one scan's GroupFit result back the way FitModel reads a single fit: summary, curve, data, averages."""
        translator = self._translator
        fit_result = result.fits[y_col]
        # Evaluated on a fine grid rather than at the scan's own points, as FitModel draws a single fit.
        x_fine = np.linspace(scan.x.min(), scan.x.max(), 300)
        curve = FitCurve(
            source_scan_uuid=scan.uuid,
            scan_name=scan.name,
            x=x_fine.tolist(),
            best_fit=np.asarray(fit_result.raw.eval(x=x_fine)).tolist(),
            components=translator._evaluate_components(fit_result, x_fine),
        )
        return GroupFitOutcome(
            source_scan_uuid=scan.uuid,
            spec=scan.spec,
            result=translator._build_result_summary(scan.spec, fit_result),
            curve=curve,
            data=FitData(x=scan.x.tolist(), y=scan.y.tolist(), err=np.sqrt(np.abs(scan.y)).tolist()),
            averages={
                name: ColumnAverageSummary(mean=average.mean, std=average.std)
                for name, average in result.averages.items()
            },
        )
