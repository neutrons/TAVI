import numpy as np
import pytest

from tavi.library.fit import CosineModel, Fit, FitPackage, ModelName, cosine


def _cosine_data(amplitude=3.0, frequency=2.5, phase=0.4, offset=7.0, n=200):
    """A clean cosine A*cos(w*(x + phi)) + d0 sampled over several periods."""
    x = np.linspace(0.0, 10.0, n)
    y = amplitude * np.cos(frequency * (x + phase)) + offset
    return x, y


def _phase_error(fitted, true, frequency):
    """Phase shifts that differ by a whole period describe the same curve."""
    return abs((frequency * (fitted - true) + np.pi) % (2 * np.pi) - np.pi)


def test_cosine_evaluates_the_stated_formula():
    """The model function is A*cos(w*(x + phi)) + d0."""
    x = np.linspace(-1.0, 1.0, 11)
    expected = 3.0 * np.cos(2.5 * (x + 0.4)) + 7.0
    assert np.allclose(cosine(x, amplitude=3.0, frequency=2.5, phase=0.4, offset=7.0), expected)


def test_cosine_guess_lands_near_the_true_parameters():
    """guess reads frequency off a periodogram, so it starts near the right period."""
    x, y = _cosine_data()
    guess = Fit(package=FitPackage.lmfit).guess(x, y, ModelName.Cosine)
    assert guess["amplitude"] == pytest.approx(3.0, abs=0.1)
    assert guess["frequency"] == pytest.approx(2.5, abs=0.1)
    assert guess["offset"] == pytest.approx(7.0, abs=0.1)
    assert _phase_error(guess["phase"], 0.4, 2.5) < 0.2


def test_cosine_fit_recovers_all_four_parameters():
    """A fit seeded from guess recovers A, w, phi and d0."""
    x, y = _cosine_data()
    result = Fit(package=FitPackage.lmfit).fit(x, y, [(ModelName.Cosine, dict(guess=True))])
    assert result.amplitude == pytest.approx(3.0, abs=1e-3)
    assert result.frequency == pytest.approx(2.5, abs=1e-3)
    assert result.offset == pytest.approx(7.0, abs=1e-3)
    assert _phase_error(result.phase, 0.4, 2.5) < 1e-3


def test_cosine_fit_recovers_parameters_on_uneven_x():
    """Unevenly spaced x needs no resampling - the frequency guess is Lomb-Scargle."""
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0.0, 10.0, 150))
    y = 3.0 * np.cos(2.5 * (x + 0.4)) + 7.0
    result = Fit(package=FitPackage.lmfit).fit(x, y, [(ModelName.Cosine, dict(guess=True))])
    assert result.frequency == pytest.approx(2.5, abs=1e-3)
    assert result.amplitude == pytest.approx(3.0, abs=1e-3)


def test_cosine_fit_with_explicit_initial_parameters():
    """All four cosine parameters can be seeded explicitly instead of guessed."""
    x, y = _cosine_data()
    result = Fit(package=FitPackage.lmfit).fit(
        x,
        y,
        [(ModelName.Cosine, dict(amplitude=2.0, frequency=2.4, phase=0.0, offset=6.0))],
    )
    assert result.frequency == pytest.approx(2.5, abs=1e-3)
    assert result.offset == pytest.approx(7.0, abs=1e-3)


def test_cosine_partial_initial_parameters_raise():
    """A half-specified cosine is an error rather than a silent default."""
    x, y = _cosine_data()
    with pytest.raises(ValueError, match="not supported"):
        Fit(package=FitPackage.lmfit).fit(x, y, [(ModelName.Cosine, dict(amplitude=2.0, frequency=2.4))])


def test_cosine_frequency_can_be_fixed():
    """set works on cosine parameters like any other model's."""
    x, y = _cosine_data()
    result = Fit(package=FitPackage.lmfit).fit(
        x,
        y,
        [(ModelName.Cosine, dict(guess=True, set={"frequency": dict(value=1.0, vary=False)}))],
    )
    assert result.frequency == 1.0


def test_cosine_composes_with_another_component():
    """A cosine carries a prefix and composes into a multi-component fit."""
    x, y = _cosine_data()
    y = y + 0.5 * x
    result = Fit(package=FitPackage.lmfit).fit(
        x,
        y,
        [
            (ModelName.Cosine, dict(guess=True, prefix="cos_")),
            (ModelName.Linear, dict(guess=True, prefix="bg_")),
        ],
    )
    assert result["cos_"].values["frequency"] == pytest.approx(2.5, abs=1e-2)
    assert result["bg_"].values["slope"] == pytest.approx(0.5, abs=1e-2)
    # No component has a center, so the fit reports no peaks.
    assert result.peaks == []


def test_cosine_guess_requires_x():
    """The frequency estimate is meaningless without the independent variable."""
    _, y = _cosine_data()
    with pytest.raises(ValueError, match="requires the independent variable"):
        CosineModel().guess(y)


def test_cosine_guess_on_flat_data_falls_back_to_one_cycle():
    """Flat data has no periodogram peak; the guess is one cycle across the range."""
    x = np.linspace(0.0, 10.0, 50)
    guess = Fit(package=FitPackage.lmfit).guess(x, np.full_like(x, 4.0), ModelName.Cosine)
    assert guess["frequency"] == pytest.approx(2 * np.pi / 10.0)
    assert guess["amplitude"] == 0.0
    assert guess["offset"] == pytest.approx(4.0)
