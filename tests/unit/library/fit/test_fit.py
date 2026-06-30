import numpy as np
import pytest

from tavi.library.fit import Fit, FitPackage, ModelName


def _gaussian_data():
    """A clean Gaussian centered at 0.5 on a flat-ish baseline."""
    x = np.linspace(-2, 3, 200)
    y = 100 * np.exp(-((x - 0.5) ** 2) / (2 * 0.2**2))
    return x, y


def test_set_fixes_parameter():
    """A parameter with vary=False stays at its set value through the fit."""
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(
        x,
        y,
        [(ModelName.Gaussian, dict(guess=True, set={"center": dict(value=0.0, vary=False)}))],
    )
    # Center was frozen away from the true 0.5, so it must remain exactly 0.0.
    assert result.peak.values["center"] == 0.0
    assert result.raw.params["center"].vary is False


def test_set_bounds_parameter():
    """Bounds from set are honored by the fit."""
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(
        x,
        y,
        [(ModelName.Gaussian, dict(guess=True, set={"sigma": dict(max=0.1)}))],
    )
    assert result.peak.values["sigma"] <= 0.1 + 1e-9
    assert result.raw.params["sigma"].max == 0.1


def test_set_uses_prefix():
    """set keys are bare names; the component prefix is added automatically."""
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(
        x,
        y,
        [(ModelName.Gaussian, dict(guess=True, prefix="g1_", set={"center": dict(value=0.0, vary=False)}))],
    )
    assert result["g1_"].values["center"] == 0.0
    assert result.raw.params["g1_center"].vary is False


def test_explicit_linear_params():
    """A linear background can be seeded with explicit slope/intercept."""
    x = np.linspace(-2, 3, 200)
    y = 2.0 * x + 5.0
    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(x, y, [(ModelName.Linear, dict(slope=1.0, intercept=0.0))])
    assert np.isclose(result.components[""].values["slope"], 2.0)
    assert np.isclose(result.components[""].values["intercept"], 5.0)


def test_set_unknown_parameter_raises():
    """Setting a parameter that does not exist is an explicit error."""
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    with pytest.raises(ValueError, match="unknown parameter"):
        fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True, set={"nope": dict(value=1.0)}))])
