"""Fit several columns of every scan in a group against one shared x column, and average others."""

import copy
from dataclasses import dataclass, field, replace
from typing import Any, Optional, Sequence

import numpy as np

from tavi.library.data.scan import UUID, Scan
from tavi.library.fit.fit import Fit, FitPackage, FitResult, ModelName

ModelDict = list[tuple[ModelName, dict[str, Any]]]

# Starting values Fit.fit accepts, across every supported model shape. Chaining copies only these from the
# previous scan's fit, since derived values (fwhm, height) and Voigt's gamma are not Fit.fit inputs.
_SEED_KEYS = ("center", "sigma", "amplitude", "slope", "intercept", "frequency", "phase", "offset")

# A seed smaller than this fraction of its parameter's data scale is snapped to 0 - see _snap_zero.
_SEED_ZERO_FRACTION = 1e-6


@dataclass
class GroupFitParams:
    """
    Which columns of every scan to fit against a shared x column, and which to average.

    Attributes:
        x_col: Independent-variable column; None uses each scan's ``tavimeta.default_axis[0]``.
        fits: y column -> composite model in the form :meth:`Fit.fit` takes. A peak and a sinusoid
            differ only in their model (e.g. ``ModelName.Gaussian`` vs ``ModelName.Cosine``).
        averages: Columns reduced to their per-scan mean and standard deviation instead of fit.
        overrides: scan uuid -> y column -> component prefix -> patch merged into that component of
            ``fits[y column]`` for that scan only, e.g. to fix a weak peak's sigma in one scan.
            ``set`` options merge per parameter; every other key replaces the group value.
        chain: Start each scan after the first from the previous scan's fitted values instead of the
            group model's start/guess, so a peak that drifts or fades is tracked from its neighbour's
            solution. ``set`` options and overrides still apply on top of the seed.

    """

    x_col: Optional[str] = None
    fits: dict[str, ModelDict] = field(default_factory=dict)
    averages: list[str] = field(default_factory=list)
    overrides: dict[UUID, dict[str, dict[str, dict[str, Any]]]] = field(default_factory=dict)
    chain: bool = True


@dataclass
class ColumnAverage:
    """Per-scan mean of a column and its standard deviation (its spread during the scan)."""

    mean: float
    std: float


@dataclass
class GroupFitResult:
    """One scan's results: its uuid, a FitResult per fitted column and a ColumnAverage per averaged column."""

    uuid: UUID
    fits: dict[str, FitResult]
    averages: dict[str, ColumnAverage]


class GroupFit:
    """Apply one set of per-column fits and averages to every scan of a group, with per-scan fit overrides."""

    def __init__(
        self,
        scans: Sequence[Scan],
        group_fit_params: GroupFitParams,
        package: FitPackage = FitPackage.lmfit,
    ) -> None:
        """Initialize with the scans to process, in order, what to do with their columns, and the fit backend."""
        self.scans = list(scans)
        self.results: list[GroupFitResult] = []
        self._fit = Fit(package=package)
        # Keep a private copy so override() never reaches back into the caller's params, and route
        # up-front overrides through override() so they are validated the same way as later ones.
        self.group_fit_params = replace(copy.deepcopy(group_fit_params), overrides={})
        for uuid, by_column in group_fit_params.overrides.items():
            for y_col, by_prefix in by_column.items():
                for prefix, patch in by_prefix.items():
                    self.override(uuid, y_col, prefix, patch)

    def fit(self) -> list[GroupFitResult]:
        """Fit and average every scan's requested columns, returning one :class:`GroupFitResult` per scan."""
        params = self.group_fit_params
        if not params.fits and not params.averages:
            raise ValueError("group_fit_params has no columns to fit or average.")

        self.results = []
        for scan in self.scans:
            averages = {}
            for name in params.averages:
                values = self._column(scan, name)
                averages[name] = ColumnAverage(mean=float(np.mean(values)), std=float(np.std(values)))
            x = self._column(scan, params.x_col or scan.tavimeta.default_axis[0])
            fits = {}
            for y_col in params.fits:
                y = self._column(scan, y_col)
                # Counting statistics on the raw column, matching how single-scan fits weight points.
                fits[y_col] = self._fit.fit(x, y, self.model_for(scan.uuid, y_col), err=np.sqrt(np.abs(y)))
            self.results.append(GroupFitResult(uuid=scan.uuid, fits=fits, averages=averages))
        return self.results

    def model_for(self, uuid: UUID, y_col: str) -> ModelDict:
        """Return the model ``y_col`` is fit with in scan ``uuid``: group model, chained seed, then that scan's patches."""
        params = self.group_fit_params
        index = [scan.uuid for scan in self.scans].index(uuid)
        # fit() appends results as it goes, so the previous scan's result exists by the time this scan is fit.
        previous = self.results[index - 1].fits.get(y_col) if params.chain and 0 < index <= len(self.results) else None
        patches = params.overrides.get(uuid, {}).get(y_col, {})
        model_dict = []
        for model_name, initial_params in params.fits[y_col]:
            prefix = initial_params.get("prefix", "")
            if previous is not None and prefix in previous.components:
                seed = {key: value for key, value in previous[prefix].values.items() if key in _SEED_KEYS}
                # A non-finite seed would poison this fit; keep the component's own start instead.
                if seed and all(np.isfinite(value) for value in seed.values()):
                    scales = self._seed_scales(self.scans[index], y_col, model_name)
                    seed = {key: _snap_zero(value, scales[key]) for key, value in seed.items()}
                    initial_params = {**{k: v for k, v in initial_params.items() if k != "guess"}, **seed}
            if prefix in patches:
                initial_params = _merge(initial_params, patches[prefix])
            model_dict.append((model_name, initial_params))
        return model_dict

    def override(self, uuid: UUID, y_col: str, prefix: str, patch: dict[str, Any]) -> None:
        """Merge ``patch`` into component ``prefix`` of ``y_col``'s model for scan ``uuid`` only; applies on the next fit()."""
        if uuid not in [scan.uuid for scan in self.scans]:
            raise KeyError(f"Scan {uuid.value!r} is not in this group.")
        prefixes = [initial_params.get("prefix", "") for _, initial_params in self.group_fit_params.fits[y_col]]
        if prefix not in prefixes:
            raise KeyError(f"No component {prefix!r} in the {y_col!r} model; available: {prefixes}.")
        by_prefix = self.group_fit_params.overrides.setdefault(uuid, {}).setdefault(y_col, {})
        by_prefix[prefix] = _merge(by_prefix.get(prefix, {}), patch)

    def _seed_scales(self, scan: Scan, y_col: str, model_name: ModelName) -> dict[str, float]:
        """Return each seedable parameter's natural size in ``scan``'s data, for a component of ``model_name``."""
        x = self._column(scan, self.group_fit_params.x_col or scan.tavimeta.default_axis[0])
        y = self._column(scan, y_col)
        x_scale = float(np.max(np.abs(x))) if x.size else 0.0
        x_span = float(np.ptp(x)) if x.size else 0.0
        y_scale = float(np.max(np.abs(y))) if y.size else 0.0
        return {
            "center": x_scale,
            "sigma": x_span,
            # A peak's amplitude is its integrated area; a cosine's is a height.
            "amplitude": y_scale if model_name is ModelName.Cosine else y_scale * x_span,
            "slope": y_scale / x_span if x_span else y_scale,
            "intercept": y_scale,
            # One cycle across the scan, the scale CosineModel.guess falls back to.
            "frequency": 2 * np.pi / x_span if x_span else 1.0,
            "phase": x_span,
            "offset": y_scale,
        }

    @staticmethod
    def _column(scan: Scan, name: str) -> np.ndarray:
        """Return a scan column as floats, accepting the leading-underscore spelling of digit-leading names."""
        data = scan.data.data
        if name not in data and name[:1].isdigit() and f"_{name}" in data:
            name = f"_{name}"
        if name not in data:
            raise KeyError(f"Column {name!r} not found in scan {scan.uuid.value}; available: {list(data)}.")
        return np.asarray(data[name], dtype=float)


def _snap_zero(value: float, scale: float) -> float:
    """
    Return ``value``, or exactly 0 when it is negligible next to ``scale``.

    A fit that lands on (numerically) zero, e.g. a center of 2e-13, would otherwise seed the next
    scan with it - and lmfit's finite-difference step is proportional to the starting value, so
    that parameter stays frozen for the rest of the chain. Starting from exactly 0 makes lmfit
    fall back to an absolute step, as the deprecated sequential fit's ``FitModel._seeded_value`` did.
    """
    return 0.0 if abs(value) < _SEED_ZERO_FRACTION * scale else value


def _merge(params: dict[str, Any], patch: dict[str, Any]) -> dict[str, Any]:
    """Return ``params`` with ``patch`` applied: ``set`` options merge per parameter, other keys replace."""
    base_set = params.get("set", {})
    merged = {**params, **copy.deepcopy(patch)}
    merged["set"] = {
        **base_set,
        **{name: {**base_set.get(name, {}), **opts} for name, opts in patch.get("set", {}).items()},
    }
    return merged
