"""Hb1a specific plugins to reduce diffraction data."""

import warnings
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
        """
        Sum a fit's peak amplitudes, combining their errors in quadrature (assuming independent peaks).

        A fit lmfit could not measure contributes (0, 0): it drops every stderr when it cannot
        invert the Hessian, and an unmeasured reflection is still exported, carrying no weight.
        """
        peaks = fit_result.peaks
        if not peaks:
            return 0.0, 0.0
        if any(peak.values.get("amplitude") is None or peak.errors.get("amplitude") is None for peak in peaks):
            return 0.0, 0.0
        amplitude = sum(peak.values["amplitude"] for peak in peaks)
        amplitude_err = np.sqrt(sum(peak.errors["amplitude"] ** 2 for peak in peaks))
        return amplitude, amplitude_err

    @staticmethod
    def remove_peak(scan_list: list[int], fit_results: list, no_peak: list[int]) -> list:
        """
        Zero the peak amplitude and error of each no_peak scan, modifying and returning fit_results.

        scan_list and fit_results are positional, as browse returns them. A scan measured
        where no peak turned out to be stays in the export carrying zero intensity, since a
        refinement cannot otherwise tell it from one that was never measured. A no_peak scan
        missing from scan_list is warned about rather than raised on.
        """
        if len(scan_list) != len(fit_results):
            raise ValueError(
                f"scan_list and fit_results must be the same length, got {len(scan_list)} and {len(fit_results)}."
            )

        missing = [scan for scan in no_peak if scan not in scan_list]
        if missing:
            warnings.warn(
                f"no_peak scans {missing} are not in scan_list, so they have no fit result to zero.",
                stacklevel=2,
            )

        zeroed = set(no_peak)
        for scan, fit_result in zip(scan_list, fit_results):
            if scan not in zeroed:
                continue
            # Only the peak components are zeroed; a background component carries no
            # amplitude, and keeping it leaves the fit still plottable against its data.
            for component in fit_result.peaks:
                component.values["amplitude"] = 0.0
                component.errors["amplitude"] = 0.0
        return fit_results

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
        magnetic: bool = False,
    ) -> None:
        """
        Write (hkl, intensity, error) entries to save_to_file in the .int format used for refinement.

        magnetic writes each line's code in an i4 field between the hkl and the intensity,
        telling the refinement which reflection a line belongs to when several share an hkl
        (default 1), and carries the declared format line with it so the header cannot
        disagree with the lines; without it no code field is written and codes is unused.
        Every line ends in the refinement weight, always 1 for these data.
        propagation_vectors declares the k vectors the codes index, after the wavelength
        line; None omits the block.
        """
        target = save_to_file if overwrite else VERITAS._next_version_path(save_to_file)
        codes = codes if codes is not None else [1] * len(export)
        with open(target, "w") as f:
            f.write(f"{title}\n")
            f.write("(3i5,i4,2f8.2,i4)\n" if magnetic else "(3i5,2f8.2,i4)\n")
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
                code_field = f"{int(code):4d}" if magnetic else ""
                # The trailing i4 is the refinement weight, which every measured
                # reflection here carries equally.
                f.write(f"{hkl_fields}{code_field}{intensity:8.2f}{err:8.2f}{1:4d}\n")

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
        wavevector: Optional[list] = None,
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
            overwrite: True (default) replaces save_to_file; False writes "_<n>" before the
                suffix, n being the next unused version on disk.
            background_results: Background run per peak, one per entry in hkls. Each peak's
                amplitude has its background subtracted and the errors added in quadrature.
                Left empty (default), amplitudes are exported as fitted.
            wavevector: The propagation vector declared in the header when background_results
                is given, each line's code field indexing it. Required to write such a file.

        """
        if background_results and len(background_results) != len(hkls):
            raise ValueError(
                f"background_results must be one fit result per peak ({len(hkls)}), got {len(background_results)}."
            )
        if save_to_file and background_results and wavevector is None:
            raise ValueError("background_results subtracts a nuclear run, so its wavevector must be given.")

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
            # A background-subtracted export is the magnetic scattering left over, so the
            # refinement needs the propagation vector its lines index declared in the header.
            # The default codes are already 1, which is that vector.
            VERITAS._write_int_file(
                save_to_file,
                overwrite,
                title,
                wavelength,
                export,
                propagation_vectors=[wavevector] if background_results else None,
                magnetic=bool(background_results),
            )
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

        A line indexes its k vector by the absolute value of its code, a negative code meaning
        H - k, so branches sharing a vector number must resolve to the same k.
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

        Each satellite is written at its parent nuclear reflection, the +q branch shifted by
        ``hkl - wavevector`` and the -q branch by ``hkl + wavevector``, so the returned list
        alternates the two. The branches need not be the same length - where one runs out, the
        remainder of the longer is appended on its own. Both branches land on the same hkl, so
        only the satellite code between each line's hkl and intensity tells them apart; the k
        vectors those codes index are declared in the header.

        Args:
            title: Title line written as the first line of the file header.
            plus: The +q branch, as ``[hkls, fit_results, res_4ds, wavevector]`` - one (h, k, l),
                fit result and (resolution matrix, r0) per peak, plus the branch's wavevector.
            minus: The -q branch, laid out as plus. It need not hold the same number of peaks.
            ax: Scan axis; "s1" for transverse scans.
            save_to_file: Output file path. No file is written if None.
            wavelength: Neutron wavelength written to the file header.
            overwrite: True (default) replaces save_to_file; False writes "_<n>" before the
                suffix, n being the next unused version on disk.
            satellite_codes: The (plus, minus) i4 values telling the branches apart. (1, 2)
                gives each its own propagation vector; (1, -1) makes the second a signed index
                into one vector, requiring both branches to share a wavevector.

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
                magnetic=True,
            )
        return export
