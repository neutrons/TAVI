"""Fit module."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

import lmfit
import numpy as np


class FitPackage(Enum):
    """Supported fitting backends."""

    lmfit = "lmfit"


class ModelName(Enum):
    """Supported peak/background model shapes."""

    Gaussian = "Gaussian"
    Lorentzian = "Lorentzian"
    Voigt = "Voigt"
    Custom = "Custom defined model"
    Linear = "linear model"

# @dataclass
# class FitData:
#     """This will be in the GUI. It will have Provenance(origin of data, normalization, limit, column names). uuid, FitResult."""

@dataclass
class ComponentResult:
    """
    Fitted parameters for a single model component, identified by its prefix.

    A composite fit is made of one or more components (e.g. ``g1_``, ``g2_``,
    ``exp_``). Each component stores the bare parameter names (with the prefix
    stripped) so callers can read e.g. ``component["center"]`` regardless of the
    prefix used in the composite model.

    Attributes:
        prefix: The prefix lmfit assigned to this component ("" if none).
        values: Bare parameter name -> fitted value (e.g. center, sigma, amplitude).
        errors: Bare parameter name -> 1-sigma uncertainty, or None if not estimated.

    """

    prefix: str
    values: dict[str, float]
    errors: dict[str, Optional[float]]

    def __getitem__(self, key: str) -> float:
        """Return the fitted value of parameter ``key`` (without the prefix)."""
        return self.values[key]


@dataclass
class FitResult:
    """
    Backend-agnostic result of a (possibly multi-component) peak fit.

    Holds one :class:`ComponentResult` per model component, keyed by prefix, so
    callers do not depend on any specific fitting library. The native backend
    result is kept in ``raw`` as an escape hatch for advanced use.

    For a single-component fit, the common scalar parameters (``center``,
    ``sigma``, ``amplitude``, ``fwhm``, ``height``) and their ``*_err``
    counterparts can be read directly as attributes; they delegate to the only
    component. With multiple components, read them via ``result[prefix]``.

    Attributes:
        components: Prefix -> :class:`ComponentResult` for every model component.
        redchi: Reduced chi-square of the fit.
        best_fit: Model evaluated at the input x (same shape as the data).
        raw: The native backend result object.
        fit_function: The (composite) model object used for the fit.

    """

    components: dict[str, ComponentResult]
    redchi: float
    best_fit: np.ndarray
    raw: Any
    fit_function: Any

    def __getitem__(self, prefix: str) -> ComponentResult:
        """Return the component fitted with the given ``prefix``."""
        return self.components[prefix]

    def __getattr__(self, name: str) -> Any:
        """Delegate scalar parameter access to the sole component when unambiguous."""
        components = self.__dict__.get("components", {})
        if len(components) == 1:
            (component,) = components.values()
            if name.endswith("_err"):
                base = name[:-4]
                if base in component.errors:
                    return component.errors[base]
            elif name in component.values:
                return component.values[name]
        raise AttributeError(name)


class Fit:
    """Provide a universal interface for tavi fit."""

    def __init__(self, package: FitPackage = FitPackage.lmfit) -> None:
        """Initialize with the fitting backend to use."""
        self.package = package

    def fit(
        self,
        x: np.ndarray,
        y: np.ndarray,
        model_dict: list[tuple[ModelName, dict[str, Any]]],
    ) -> FitResult:
        """
        Fit a composite model built from one or more named sub-models.

        Each sub-model may carry a ``prefix`` so its parameters are namespaced in
        the composite fit (e.g. ``g1_center``, ``g2_center``, ``exp_decay``).
        Currently only the lmfit backend is supported.

        Args:
            x: Independent variable values.
            y: Measured values to fit.
            model_dict: Sequence of ``(model_name, initial_params)`` pairs.
                ``initial_params`` is a dict that may hold:
                  - ``prefix``: parameter namespace for this component (default "").
                  - ``guess``: if truthy, let the model guess initial parameters
                    from the data instead of using explicit values.
                  - ``center`` / ``sigma`` / ``amplitude``: initial values used
                    when ``guess`` is not set.

        Returns:
            A :class:`FitResult` with one :class:`ComponentResult` per component,
            keyed by prefix.

        """
        if self.package != FitPackage.lmfit:
            raise ValueError(f"Package {self.package} not supported yet.")

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        fit_function = None
        params = lmfit.Parameters()
        prefixes = []
        for model_name, initial_params in model_dict:
            prefix = initial_params.get("prefix", "")
            prefixes.append(prefix)
            model = self._build_model(model_name, prefix)
            if initial_params.get("guess", False):
                params.update(model.guess(y, x=x))
            else:
                params.update(
                    model.make_params(
                        center=initial_params["center"],
                        sigma=initial_params["sigma"],
                        amplitude=initial_params["amplitude"],
                    )
                )
            fit_function = model if fit_function is None else fit_function + model

        if fit_function is None:
            raise ValueError("model_dict is empty; provide at least one model.")

        result = fit_function.fit(y, params, x=x)
        return self._build_fit_result(result, prefixes, fit_function)

    @staticmethod
    def _build_model(model_name: ModelName, prefix: str) -> lmfit.Model:
        """Construct an lmfit model for the given model name and prefix."""
        match model_name:
            case ModelName.Gaussian:
                return lmfit.models.GaussianModel(prefix=prefix)
            case ModelName.Lorentzian:
                return lmfit.models.LorentzianModel(prefix=prefix)
            case ModelName.Voigt:
                return lmfit.models.VoigtModel(prefix=prefix)
            case ModelName.Linear:
                return lmfit.models.LinearModel(prefix=prefix)
            case ModelName.Custom:
                raise ValueError("Need to think about how to parse a custom defined function in lmfit.")
            case _:
                raise ValueError("Model not recognized.")

    @staticmethod
    def _build_fit_result(result: Any, prefixes: list[str], fit_function: Any) -> FitResult:
        """Group fitted parameters by prefix into a FitResult."""
        components: dict[str, ComponentResult] = {}
        for prefix in prefixes:
            values: dict[str, float] = {}
            errors: dict[str, Optional[float]] = {}
            for name, par in result.params.items():
                # Assign each parameter to its longest matching prefix so an
                # empty-prefix component does not swallow prefixed parameters.
                owners = [p for p in prefixes if name.startswith(p)]
                if not owners or max(owners, key=len) != prefix:
                    continue
                base = name[len(prefix) :]
                values[base] = par.value
                errors[base] = par.stderr
            components[prefix] = ComponentResult(prefix=prefix, values=values, errors=errors)

        return FitResult(
            components=components,
            redchi=result.redchi,
            best_fit=result.best_fit,
            raw=result,
            fit_function=fit_function,
        )
