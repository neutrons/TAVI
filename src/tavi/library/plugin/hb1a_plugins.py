"""Hb1a specific plugins to reduce diffraction data."""

from pathlib import Path
from typing import Optional

import numpy as np

from tavi.library.fit import FitResult


class VERITAS:
    """Plugins."""

    @staticmethod
    def _next_version_path(file_path: str) -> str:
        """Return file_path with "_<n>" inserted before the suffix, n being the next unused version."""
        path = Path(file_path)
        version = 1
        candidate = path.with_name(f"{path.stem}_{version}{path.suffix}")
        while candidate.exists():
            version += 1
            candidate = path.with_name(f"{path.stem}_{version}{path.suffix}")
        return str(candidate)

    @staticmethod
    def _summed_amplitude(fit_result: FitResult) -> tuple[float, float]:
        """Sum a fit's peak amplitudes, combining their errors in quadrature (assuming independent peaks)."""
        peaks = fit_result.peaks
        amplitude = sum(peak.values["amplitude"] for peak in peaks)
        amplitude_err = np.sqrt(sum(peak.errors["amplitude"] ** 2 for peak in peaks))
        return amplitude, amplitude_err

    @staticmethod
    def _lorentz_factor(res_4d: list, ax: str) -> float:
        """Resolution-derived Lorentz factor of one peak, for a scan along ax."""
        mat, r0 = res_4d[0][0], res_4d[0][1]
        det = np.abs(mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0])
        if ax == "s1":
            # use mat[1, 1] for transverse scans, R0 factor is optional
            return r0 * np.sqrt(det) / np.sqrt(mat[1, 1]) / np.sqrt(2 * np.pi)
        if ax == "th2th":
            return r0 * np.sqrt(det) / np.sqrt(mat[0, 0]) / np.sqrt(2 * np.pi)
        raise ValueError("axis not defined. Only supporting export s1 and the2th scans currently.")

    @staticmethod
    def _write_int_file(
        save_to_file: str,
        overwrite: bool,
        title: str,
        wavelength: float,
        export: list,
        codes: Optional[list] = None,
        propagation_vectors: Optional[list] = None,
        code_before_intensity: bool = False,
    ) -> None:
        """
        Write (hkl, intensity, error) entries to save_to_file in the .int format used for refinement.

        codes is the integer written in each line's i4 field; it tells the refinement which
        reflection a line belongs to when several share an hkl. Defaults to 1 for every
        line, as a commensurate export has nothing to distinguish.

        code_before_intensity puts that field between the hkl and the intensity rather than
        after the error. It moves the declared format line with it, since the two describe
        the same columns and a file whose header disagreed with its lines would be misread.

        propagation_vectors declares the k vectors those codes index, written after the
        wavelength line as a count followed by one "h k l" line each. Left as None, the
        block is omitted entirely, as a commensurate export has no k vector to declare.
        """
        target = save_to_file if overwrite else VERITAS._next_version_path(save_to_file)
        codes = codes if codes is not None else [1] * len(export)
        with open(target, "w") as f:
            f.write(f"{title}\n")
            f.write("(3i5,i4,2f8.2,3f8.2)\n" if code_before_intensity else "(3i5,2f8.2,i4,3f8.2)\n")
            f.write(f"{wavelength}  0   0\n")

            if propagation_vectors is not None:
                f.write(f"{len(propagation_vectors)}\n")
                for vector in propagation_vectors:
                    # .10g keeps an incommensurate component exact enough to refine
                    # against while still writing a whole number as "1" rather than "1.0".
                    # Negating a k vector turns its zero components into -0.0, so they are
                    # folded back to 0.0 rather than written out as "-0".
                    f.write(" ".join(f"{component if component else 0.0:.10g}" for component in vector) + "\n")

            for ((h, k, l), intensity, err), code in zip(export, codes):
                hkl_fields = f"{int(round(h)):5d}{int(round(k)):5d}{int(round(l)):5d}"
                if code_before_intensity:
                    f.write(f"{hkl_fields}{int(code):4d}{intensity:8.2f}{err:8.2f}\n")
                else:
                    f.write(f"{hkl_fields}{intensity:8.2f}{err:8.2f}{int(code):4d}\n")

    @staticmethod
    def export_intensity(
        title: str,
        hkls: list,
        fit_results: list,
        res_4ds: list,
        ax: str,
        save_to_file: Optional[str],
        wavelength: float = 2.37815,
        overwrite: bool = True,
        background_results: Optional[list] = [],
    ) -> list:
        """
        Export the intensity data to a .int file for refinement.

        Args:
            title: Title line written as the first line of the file header.
            hkls: List of (h, k, l) per peak.
            fit_results: Fit result per peak, providing amplitude and amplitude_err.
            res_4ds: (resolution matrix, r0) per peak.
            ax: Scan axis; "s1" for transverse scans.
            save_to_file: Output file path. No file is written if None.
            wavelength: Neutron wavelength written to the file header.
            overwrite: If True (default), write to save_to_file, replacing it if it
                exists. If False, write to a new file with "_<n>" appended before the
                suffix, where n is the next unused version number on disk.
            background_results: Fit result of the background run measured for each peak,
                one per entry in hkls. When given, each peak's amplitude has its
                background's amplitude subtracted, and the two errors are combined in
                quadrature. Left empty (default), the amplitudes are exported as fitted.

        """
        if background_results and len(background_results) != len(hkls):
            raise ValueError(
                f"background_results must be one fit result per peak ({len(hkls)}), got {len(background_results)}."
            )

        # zip stops at the shortest input, so an absent background list has to be padded
        # rather than left empty - otherwise no peak would be exported at all.
        backgrounds = background_results if background_results else [None] * len(hkls)
        export: list = []
        for hkl, fit_result, res_4d, background_result in zip(hkls, fit_results, res_4ds, backgrounds):
            lorentz_factor = VERITAS._lorentz_factor(res_4d, ax)
            # Sum the amplitudes of all peak components; combine their
            # amplitude errors in quadrature (assuming independent peaks).
            amplitude, amplitude_err = VERITAS._summed_amplitude(fit_result)

            # The background is an independent measurement, so subtracting it leaves the
            # difference less precise than either run: the errors add in quadrature even
            # though the amplitudes subtract.
            if background_result is not None:
                bkg_amplitude, bkg_amplitude_err = VERITAS._summed_amplitude(background_result)
                amplitude = amplitude - bkg_amplitude
                amplitude_err = np.sqrt(amplitude_err**2 + bkg_amplitude_err**2)

            intensity = amplitude / lorentz_factor
            err = amplitude_err / lorentz_factor
            export.append((hkl, intensity, err))
        if save_to_file:
            VERITAS._write_int_file(save_to_file, overwrite, title, wavelength, export)
        return export

    @staticmethod
    def _satellite_entry(hkl: list, fit_result: FitResult, res_4d: list, shift: np.ndarray, ax: str) -> tuple:
        """Reduce one satellite peak to its (parent hkl, intensity, error) export entry."""
        lorentz_factor = VERITAS._lorentz_factor(res_4d, ax)
        # Sum the amplitudes of all peak components; combine their
        # amplitude errors in quadrature (assuming independent peaks).
        amplitude, amplitude_err = VERITAS._summed_amplitude(fit_result)
        return (
            np.asarray(hkl, dtype=float) + shift,
            amplitude / lorentz_factor,
            amplitude_err / lorentz_factor,
        )

    @staticmethod
    def _propagation_vectors(coded_vectors: tuple) -> list:
        """
        Order the (code, k vector) pairs into the k1, k2, ... list the file header declares.

        A satellite line indexes its k vector by the absolute value of its code, and a
        negative code means the line sits at H - k of that vector, so two branches sharing
        a vector number must resolve to the same k once the sign is applied.
        """
        vectors: dict[int, np.ndarray] = {}
        for code, vector in coded_vectors:
            if code == 0:
                raise ValueError("satellite codes index propagation vectors from 1, so 0 is not a valid code.")
            resolved = vector if code > 0 else -vector
            existing = vectors.setdefault(abs(code), resolved)
            if not np.allclose(existing, resolved):
                raise ValueError(
                    f"propagation vector {abs(code)} is given two different values, {existing} and {resolved}. "
                    f"Use distinct satellite codes when the branches do not share a wavevector."
                )
        if sorted(vectors) != list(range(1, len(vectors) + 1)):
            raise ValueError(f"satellite codes must index propagation vectors 1..n, got {sorted(vectors)}.")
        return [vectors[index] for index in sorted(vectors)]

    @staticmethod
    def export_intensity_incom(
        title: str,
        plus: list,
        minus: list,
        ax: str,
        save_to_file: Optional[str],
        wavelength: float = 2.37815,
        overwrite: bool = True,
        satellite_codes: tuple[int, int] = (1, 2),
    ) -> list:
        """
        Export the satellite intensities of an incommensurate structure to a .int file for refinement.

        Each satellite is written at the parent nuclear reflection it belongs to, so the
        +q branch is shifted by ``hkl - wavevector`` and the -q branch by ``hkl + wavevector``.
        The two branches are walked together, so the returned list alternates a +q entry and
        a -q entry. The branches need not be the same length - where one runs out, the rest
        of the longer branch is appended on its own. Every entry carries its own parent hkl
        and satellite code, so a branch missing a satellite only changes the order entries
        are written in, never what any of them says.

        Both branches of a parent reflection land on the same hkl, and only the satellite
        code written in each line's i4 field, between the hkl and the intensity, tells them
        apart. The k vectors those codes index are declared in the header, after the
        wavelength line, as a count followed by one "h k l" line each.

        Args:
            title: Title line written as the first line of the file header.
            plus: The +q branch, as ``[hkls, fit_results, res_4ds, wavevector]`` - a list of
                (h, k, l) per peak, a fit result per peak providing amplitude and
                amplitude_err, a (resolution matrix, r0) per peak, and the propagation
                vector shared by the branch.
            minus: The -q branch, in the same layout as plus. It need not hold the same
                number of peaks.
            ax: Scan axis; "s1" for transverse scans.
            save_to_file: Output file path. No file is written if None.
            wavelength: Neutron wavelength written to the file header.
            overwrite: If True (default), write to save_to_file, replacing it if it
                exists. If False, write to a new file with "_<n>" appended before the
                suffix, where n is the next unused version number on disk.
            satellite_codes: The (plus, minus) integers written in the i4 field to tell the
                two branches apart. Defaults to (1, 2), declaring each branch as its own
                propagation vector; pass e.g. (1, -1) if the refinement instead expects a
                signed index into a single propagation vector, which then requires both
                branches to have been measured about the same wavevector.

        """
        hkls_plus, fit_results_plus, res_4ds_plus, wavevector_plus = plus
        hkls_minus, fit_results_minus, res_4ds_minus, wavevector_minus = minus

        if not (len(hkls_plus) == len(fit_results_plus) == len(res_4ds_plus)):
            raise ValueError("plus must hold one fit result and one resolution per hkl.")
        if not (len(hkls_minus) == len(fit_results_minus) == len(res_4ds_minus)):
            raise ValueError("minus must hold one fit result and one resolution per hkl.")

        code_plus, code_minus = satellite_codes
        # A branch's satellites sit at H + k, so shifting them back onto H shifts by -k.
        k_plus = np.asarray(wavevector_plus, dtype=float)
        k_minus = -np.asarray(wavevector_minus, dtype=float)
        branches = (
            (hkls_plus, fit_results_plus, res_4ds_plus, -k_plus, code_plus),
            (hkls_minus, fit_results_minus, res_4ds_minus, -k_minus, code_minus),
        )
        propagation_vectors = VERITAS._propagation_vectors(((code_plus, k_plus), (code_minus, k_minus)))

        export: list = []
        codes: list = []
        for index in range(max(len(hkls_plus), len(hkls_minus))):
            for hkls, fit_results, res_4ds, shift, code in branches:
                # A branch that has run out simply contributes nothing further; the other
                # one keeps going rather than its remaining satellites being dropped.
                if index >= len(hkls):
                    continue
                export.append(VERITAS._satellite_entry(hkls[index], fit_results[index], res_4ds[index], shift, ax))
                codes.append(code)
        if save_to_file:
            VERITAS._write_int_file(
                save_to_file,
                overwrite,
                title,
                wavelength,
                export,
                codes,
                propagation_vectors,
                code_before_intensity=True,
            )
        return export
