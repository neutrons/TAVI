"""Tests for the HB1A reduction plugins."""

import numpy as np
import pytest

from tavi.library.fit import ComponentResult, FitResult
from tavi.library.plugin.hb1a_plugins import VERITAS


def make_fit_result(amplitude: float, amplitude_err: float) -> FitResult:
    """Single-peak fit result carrying only the amplitude the export reads."""
    component = ComponentResult(
        prefix="p1_",
        values={"center": 0.0, "amplitude": amplitude},
        errors={"center": 0.0, "amplitude": amplitude_err},
    )
    return FitResult(
        components={"p1_": component},
        reduced_chi_squared=1.0,
        best_fit=np.zeros(1),
        raw=None,
        fit_function=None,
    )


def make_res_4d(mat: np.ndarray = None, r0: float = 1.0) -> list:
    """Resolution entry in the ``res_4d[0] == (matrix, r0)`` shape the export expects."""
    if mat is None:
        mat = np.diag([4.0, 1.0, 1.0, 1.0])
    return [(mat, r0)]


def make_branch(hkls: list, amplitudes: list, wavevector: tuple) -> list:
    """One satellite branch as ``[hkls, fit_results, res_4ds, wavevector]``."""
    fit_results = [make_fit_result(amplitude, 1.0) for amplitude in amplitudes]
    res_4ds = [make_res_4d() for _ in hkls]
    return [hkls, fit_results, res_4ds, wavevector]


def test_export_intensity_incom_shifts_each_branch_onto_parent_reflection():
    """+q satellites are shifted by -wavevector and -q satellites by +wavevector."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13), (2.0, 0.0, 0.13)], [10.0, 20.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13), (2.0, 0.0, -0.13)], [11.0, 21.0], wavevector)

    export = VERITAS.export_intensity_incom("title", plus, minus, "s1", None)

    hkls = [tuple(np.round(hkl, 6)) for hkl, _, _ in export]
    assert hkls == [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 0.0, 0.0)]


def test_export_intensity_incom_interleaves_the_two_branches():
    """Each parent reflection contributes its +q entry then its -q entry."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13), (2.0, 0.0, 0.13)], [10.0, 20.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13), (2.0, 0.0, -0.13)], [11.0, 21.0], wavevector)

    export = VERITAS.export_intensity_incom("title", plus, minus, "s1", None)

    # lorentz_factor is 1 * sqrt(4) / sqrt(1) / sqrt(2 pi) for this resolution matrix
    lorentz_factor = 2.0 / np.sqrt(2 * np.pi)
    intensities = [intensity for _, intensity, _ in export]
    assert intensities == pytest.approx([10.0, 11.0, 20.0, 21.0] / lorentz_factor)


def test_export_intensity_incom_matches_export_intensity_per_peak():
    """A satellite's intensity and error are reduced exactly as a commensurate peak's are."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [10.0], wavevector)

    incom = VERITAS.export_intensity_incom("title", plus, minus, "s1", None)
    commensurate = VERITAS.export_intensity("title", plus[0], plus[1], plus[2], "s1", None)

    for _, intensity, err in incom:
        assert intensity == pytest.approx(commensurate[0][1])
        assert err == pytest.approx(commensurate[0][2])


def test_export_intensity_incom_appends_the_longer_branchs_remainder():
    """A parent reflection measured in only one branch still gets exported, after the pairs."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13), (2.0, 0.0, 0.13), (3.0, 0.0, 0.13)], [10.0, 20.0, 30.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)

    export = VERITAS.export_intensity_incom("title", plus, minus, "s1", None)

    lorentz_factor = 2.0 / np.sqrt(2 * np.pi)
    hkls = [tuple(np.round(hkl, 6)) for hkl, _, _ in export]
    assert hkls == [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    intensities = [intensity for _, intensity, _ in export]
    assert intensities == pytest.approx([10.0, 11.0, 20.0, 30.0] / lorentz_factor)


def test_export_intensity_incom_codes_an_unpaired_satellite_by_its_own_branch(tmp_path):
    """An entry with no counterpart still carries the code of the branch it came from."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13), (2.0, 0.0, -0.13)], [11.0, 21.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target))

    lines = target.read_text().splitlines()[6:]
    assert [(line[:15], line[15:19]) for line in lines] == [
        ("    1    0    0", "   1"),
        ("    1    0    0", "   2"),
        ("    2    0    0", "   2"),
    ]


def test_export_intensity_incom_exports_an_empty_branch():
    """One branch being absent entirely leaves the other exported as measured."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13), (2.0, 0.0, 0.13)], [10.0, 20.0], wavevector)
    minus = make_branch([], [], wavevector)

    export = VERITAS.export_intensity_incom("title", plus, minus, "s1", None)

    hkls = [tuple(np.round(hkl, 6)) for hkl, _, _ in export]
    assert hkls == [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]


def test_export_intensity_incom_rejects_a_branch_missing_a_fit_result():
    """A branch's own lists must stay aligned, unlike the two branches with each other."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13), (2.0, 0.0, 0.13)], [10.0, 20.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    plus[1] = plus[1][:1]

    with pytest.raises(ValueError, match="one fit result and one resolution per hkl"):
        VERITAS.export_intensity_incom("title", plus, minus, "s1", None)


def test_export_intensity_incom_writes_parent_reflections_to_file(tmp_path):
    """The written .int file carries the shifted, integer parent hkls."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target))

    lines = target.read_text().splitlines()
    assert lines[0] == "title"
    assert len(lines) == 8  # three header lines, the k-vector block, then one line per satellite
    assert lines[6].startswith("    1    0    0")
    assert lines[7].startswith("    1    0    0")


def test_export_intensity_incom_codes_the_two_branches_apart(tmp_path):
    """Both branches share an hkl, so the i4 field carries the satellite code."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target))

    lines = target.read_text().splitlines()
    assert lines[6][15:19] == "   1"
    assert lines[7][15:19] == "   2"


def test_export_intensity_incom_writes_the_columns_its_header_declares(tmp_path):
    """The code sits between the hkl and the intensity, as the (3i5,i4,2f8.2,3f8.2) line says."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target))

    lines = target.read_text().splitlines()
    lorentz_factor = 2.0 / np.sqrt(2 * np.pi)
    assert lines[1] == "(3i5,i4,2f8.2,3f8.2)"
    assert lines[6] == f"    1    0    0   1{10.0 / lorentz_factor:8.2f}{1.0 / lorentz_factor:8.2f}"
    assert lines[7] == f"    1    0    0   2{11.0 / lorentz_factor:8.2f}{1.0 / lorentz_factor:8.2f}"


def test_export_intensity_incom_declares_both_propagation_vectors(tmp_path):
    """Codes (1, 2) give each branch its own k vector, +q for plus and -q for minus."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target))

    lines = target.read_text().splitlines()
    assert lines[3] == "2"
    assert lines[4] == "0 0 0.13"
    assert lines[5] == "0 0 -0.13"


def test_export_intensity_incom_declares_one_vector_for_signed_codes(tmp_path):
    """Codes (1, -1) make the -q branch a signed index into the single declared k vector."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target), satellite_codes=(1, -1))

    lines = target.read_text().splitlines()
    assert lines[3] == "1"
    assert lines[4] == "0 0 0.13"
    assert lines[5][15:19] == "   1"
    assert lines[6][15:19] == "  -1"


def test_export_intensity_incom_writes_an_integer_vector_without_a_decimal_point(tmp_path):
    """A whole-number k component is written as "1", matching the .int files written by hand."""
    wavevector = (0.0, 0.0, 1.0)
    plus = make_branch([(1.0, 0.0, 1.0)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -1.0)], [11.0], wavevector)
    target = tmp_path / "incom.int"

    VERITAS.export_intensity_incom("title", plus, minus, "s1", str(target), satellite_codes=(1, -1))

    assert target.read_text().splitlines()[4] == "0 0 1"


def test_export_intensity_incom_rejects_signed_codes_for_mismatched_wavevectors(tmp_path):
    """Sharing a propagation vector number requires the branches to share a wavevector."""
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], (0.0, 0.0, 0.13))
    minus = make_branch([(1.0, 0.0, -0.17)], [11.0], (0.0, 0.0, 0.17))

    with pytest.raises(ValueError, match="two different values"):
        VERITAS.export_intensity_incom("title", plus, minus, "s1", None, satellite_codes=(1, -1))


def test_export_intensity_keeps_its_own_format_line_and_column_order(tmp_path):
    """The commensurate export is unchanged by the incommensurate one: code last, no k block."""
    wavevector = (0.0, 0.0, 0.0)
    hkls, fit_results, res_4ds, _ = make_branch([(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)], [10.0, 20.0], wavevector)
    target = tmp_path / "com.int"

    VERITAS.export_intensity("title", hkls, fit_results, res_4ds, "s1", str(target))

    lines = target.read_text().splitlines()
    lorentz_factor = 2.0 / np.sqrt(2 * np.pi)
    assert lines[1] == "(3i5,2f8.2,i4,3f8.2)"
    assert len(lines) == 5  # three header lines, then one line per peak - no k-vector block
    assert lines[3] == f"    1    0    0{10.0 / lorentz_factor:8.2f}{1.0 / lorentz_factor:8.2f}   1"
    assert lines[4] == f"    2    0    0{20.0 / lorentz_factor:8.2f}{1.0 / lorentz_factor:8.2f}   1"


def test_export_intensity_incom_rejects_unknown_axis():
    """Only s1 and th2th scans have a defined Lorentz factor."""
    wavevector = (0.0, 0.0, 0.13)
    plus = make_branch([(1.0, 0.0, 0.13)], [10.0], wavevector)
    minus = make_branch([(1.0, 0.0, -0.13)], [11.0], wavevector)

    with pytest.raises(ValueError, match="axis not defined"):
        VERITAS.export_intensity_incom("title", plus, minus, "a3", None)
