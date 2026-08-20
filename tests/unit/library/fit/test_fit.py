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


# ---------------------------------------------------------------------------
# err -> weights
# ---------------------------------------------------------------------------


def test_err_downweights_a_noisy_outlier():
    """A single wildly-off point with a huge error bar must not drag the fit toward it."""
    x, y = _gaussian_data()
    y = y.copy()
    outlier_index = 10
    y[outlier_index] += 1000.0
    err = np.ones_like(x)
    err[outlier_index] = 1e6

    fit = Fit(package=FitPackage.lmfit)
    weighted = fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True))], err=err)
    unweighted = fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True))])

    true_amplitude = 100 * 0.2 * np.sqrt(2 * np.pi)
    assert abs(weighted.peak.values["amplitude"] - true_amplitude) < abs(
        unweighted.peak.values["amplitude"] - true_amplitude
    )


def test_err_zero_gets_zero_weight_not_infinite():
    """A zero-error point must not blow up into an infinite weight."""
    x, y = _gaussian_data()
    err = np.ones_like(x)
    err[0] = 0.0

    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True))], err=err)

    assert np.isfinite(result.peak.values["amplitude"])


def test_err_none_fits_unweighted_as_before():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    result = fit.fit(x, y, [(ModelName.Gaussian, dict(guess=True))], err=None)
    assert result.raw.weights is None


def test_guess_finds_center_near_true_peak():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Gaussian)
    assert guess["center"] == pytest.approx(0.5, abs=0.05)


def test_guess_amplitude_is_positive_for_a_positive_peak():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Gaussian)
    assert guess["amplitude"] > 0


def test_guess_includes_fwhm():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Gaussian)
    assert guess["fwhm"] > 0


def test_guess_strips_prefix_from_keys():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Gaussian, prefix="peak_")
    assert "center" in guess
    assert "peak_center" not in guess


def test_guess_works_for_lorentzian():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Lorentzian)
    assert guess["center"] == pytest.approx(0.5, abs=0.05)


def test_guess_works_for_voigt():
    x, y = _gaussian_data()
    fit = Fit(package=FitPackage.lmfit)
    guess = fit.guess(x, y, ModelName.Voigt)
    assert guess["center"] == pytest.approx(0.5, abs=0.05)


def test_guess_rejects_non_lmfit_package():
    x, y = _gaussian_data()
    fit = Fit(package="not-lmfit")
    with pytest.raises(ValueError, match="not supported"):
        fit.guess(x, y, ModelName.Gaussian)
