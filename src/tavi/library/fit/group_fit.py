"""Sequential fits over a group of scans, summarized against raw scan columns."""

from dataclasses import dataclass
from typing import Any, Optional, Sequence

import numpy as np

from tavi.library.data.scan import UUID, Scan
from tavi.library.fit.fit import ComponentResult, Fit, FitPackage, FitResult, ModelName

# Initial-value keys Fit.fit understands, across every supported model shape. Chaining copies only
# these from a previous fit, since derived values (fwhm, height) are not initial parameters.
_SEED_KEYS = ("center", "sigma", "amplitude", "slope", "intercept", "frequency", "phase", "offset")

_CHI_SQUARED = "reduced_chi_squared"


@dataclass
class GroupFitResult:
    """One scan's fit within a group: the scan's uuid and its FitResult."""

    uuid: UUID
    fit_result: FitResult


@dataclass
class GroupFitSummary:
    """
    One fitted parameter against one raw scan column, one point per scan.

    Attributes:
        column: Raw scan column the x axis is drawn from (e.g. ``temp``).
        parameter: Fitted parameter on the y axis (e.g. ``center``), or ``reduced_chi_squared``.
        uuids: The scan each point came from, in fit order.
        x: Per-scan mean of ``column``.
        x_err: Per-scan standard deviation of ``column`` (its spread during the scan).
        y: Fitted value of ``parameter``.
        y_err: 1-sigma uncertainty of ``parameter``; NaN where lmfit could not estimate one.

    """

    column: str
    parameter: str
    uuids: list[UUID]
    x: np.ndarray
    x_err: np.ndarray
    y: np.ndarray
    y_err: np.ndarray


class GroupFit:
    """Fit one model to each scan of a group and summarize fitted parameters against raw scan columns."""

    def __init__(self, scans: Sequence[Scan], package: FitPackage = FitPackage.lmfit) -> None:
        """Initialize with the scans to fit, in the order they should be fit, and the fitting backend."""
        self.scans = list(scans)
        self.results: list[GroupFitResult] = []
        self._fit = Fit(package=package)

    def fit(
        self,
        model_dict: list[tuple[ModelName, dict[str, Any]]],
        x_col: Optional[str] = None,
        y_col: Optional[str] = None,
        norm_col: Optional[str] = None,
        x_range: Optional[tuple[float, float]] = None,
        chain: bool = False,
    ) -> list[GroupFitResult]:
        """
        Fit ``model_dict`` to every scan in order and return one :class:`GroupFitResult` per scan.

        Args:
            model_dict: Composite model, in the form :meth:`Fit.fit` takes.
            x_col: Independent-variable column; defaults to each scan's ``tavimeta.default_axis[0]``.
            y_col: Measured column; defaults to each scan's ``tavimeta.default_axis[1]``.
            norm_col: Optional column (e.g. ``mcu``, ``time``) that y and its error are divided by.
            x_range: Optional ``(min, max)``; only points with x inside it (inclusive) are fit.
            chain: If True, seed each scan after the first with the previous scan's fitted values
                instead of the initial values/guess in ``model_dict``. Helps a sweep where a peak
                drifts or fades, so a weak scan starts from its neighbour's solution. ``set``
                options (fixed, bounded or tied parameters) still apply on top of the seed.

        """
        results: list[GroupFitResult] = []
        previous: Optional[FitResult] = None
        for scan in self.scans:
            x, y, err = self._scan_xy(scan, x_col, y_col, norm_col, x_range)
            models = self._seeded(model_dict, previous) if chain and previous is not None else model_dict
            fit_result = self._fit.fit(x, y, models, err=err)
            results.append(GroupFitResult(uuid=scan.uuid, fit_result=fit_result))
            previous = fit_result
        self.results = results
        return results

    def summarize(self, column: str, parameter: str, prefix: Optional[str] = None) -> GroupFitSummary:
        """
        Return ``parameter`` of every fit against the per-scan mean of raw ``column``.

        ``prefix`` picks the component to read; when None, the fit's single peak component is
        used (see :attr:`FitResult.peak`). ``parameter="reduced_chi_squared"`` reads fit quality.
        """
        if not self.results:
            raise ValueError("No fit results to summarize; call fit() first.")

        scans_by_uuid = {scan.uuid: scan for scan in self.scans}
        x, x_err, y, y_err = [], [], [], []
        for result in self.results:
            values = self._column(scans_by_uuid[result.uuid], column)
            x.append(np.mean(values))
            x_err.append(np.std(values))
            value, error = self._parameter(result.fit_result, parameter, prefix)
            y.append(value)
            y_err.append(error)

        return GroupFitSummary(
            column=column,
            parameter=parameter,
            uuids=[result.uuid for result in self.results],
            x=np.asarray(x, dtype=float),
            x_err=np.asarray(x_err, dtype=float),
            y=np.asarray(y, dtype=float),
            y_err=np.asarray(y_err, dtype=float),
        )

    def table(
        self, columns: Sequence[str], parameters: Sequence[str], prefix: Optional[str] = None
    ) -> dict[str, np.ndarray]:
        """Return a flat ``{name: array}`` table (uuid, columns, parameters and their ``*_err``), one row per scan."""
        if not columns and not parameters:
            raise ValueError("Ask for at least one column or parameter.")
        if not self.results:
            raise ValueError("No fit results to tabulate; call fit() first.")

        scans_by_uuid = {scan.uuid: scan for scan in self.scans}
        table: dict[str, np.ndarray] = {"uuid": np.asarray([result.uuid.value for result in self.results])}
        for column in columns:
            values = [self._column(scans_by_uuid[result.uuid], column) for result in self.results]
            table[column] = np.asarray([np.mean(v) for v in values], dtype=float)
            table[f"{column}_err"] = np.asarray([np.std(v) for v in values], dtype=float)
        for parameter in parameters:
            pairs = [self._parameter(result.fit_result, parameter, prefix) for result in self.results]
            table[parameter] = np.asarray([value for value, _ in pairs], dtype=float)
            table[f"{parameter}_err"] = np.asarray([error for _, error in pairs], dtype=float)
        return table

    @staticmethod
    def _column(scan: Scan, name: str) -> np.ndarray:
        """Return a scan column as floats, accepting the leading-underscore spelling of digit-leading names."""
        data = scan.data.data
        if name not in data and name[:1].isdigit() and f"_{name}" in data:
            name = f"_{name}"
        if name not in data:
            raise KeyError(f"Column {name!r} not found in scan {scan.uuid.value}; available: {list(data)}.")
        return np.asarray(data[name], dtype=float)

    def _scan_xy(
        self,
        scan: Scan,
        x_col: Optional[str],
        y_col: Optional[str],
        norm_col: Optional[str],
        x_range: Optional[tuple[float, float]],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Pull one scan's x, y and counting error, normalized and range-masked as requested."""
        default_x, default_y = scan.tavimeta.default_axis
        x = self._column(scan, x_col or default_x)
        y = self._column(scan, y_col or default_y)
        # Counting statistics on the raw counts, before normalization, so the error scales with y.
        err = np.sqrt(np.abs(y))
        if norm_col:
            norm = self._column(scan, norm_col)
            y = y / norm
            err = err / norm
        if x_range is not None:
            mask = (x >= min(x_range)) & (x <= max(x_range))
            x, y, err = x[mask], y[mask], err[mask]
        return x, y, err

    @staticmethod
    def _seeded(
        model_dict: list[tuple[ModelName, dict[str, Any]]], previous: FitResult
    ) -> list[tuple[ModelName, dict[str, Any]]]:
        """Return ``model_dict`` with each component's initial values replaced by ``previous``'s fitted values."""
        seeded = []
        for model_name, initial_params in model_dict:
            component = previous.components.get(initial_params.get("prefix", ""))
            seed = (
                {key: value for key, value in component.values.items() if key in _SEED_KEYS}
                if component is not None
                else {}
            )
            # A non-finite seed would poison the next fit; fall back to this component's own start.
            if not seed or not all(np.isfinite(value) for value in seed.values()):
                seeded.append((model_name, initial_params))
                continue
            params = {key: value for key, value in initial_params.items() if key != "guess"}
            params.update(seed)
            seeded.append((model_name, params))
        return seeded

    @staticmethod
    def _parameter(fit_result: FitResult, parameter: str, prefix: Optional[str]) -> tuple[float, float]:
        """Return ``(value, 1-sigma)`` of ``parameter`` from one fit, NaN error when none was estimated."""
        if parameter == _CHI_SQUARED:
            return fit_result.reduced_chi_squared, np.nan
        component: ComponentResult = fit_result.peak if prefix is None else fit_result[prefix]
        if parameter not in component.values:
            raise KeyError(
                f"Parameter {parameter!r} not in component {component.prefix!r}; available: {list(component.values)}."
            )
        error = component.errors.get(parameter)
        return component.values[parameter], np.nan if error is None else error
