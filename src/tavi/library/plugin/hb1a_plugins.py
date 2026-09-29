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
            mat, r0 = res_4d[0][0], res_4d[0][1]
            # ====================================================
            # resolution calculated here
            det = np.abs(mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0])
            if ax == "s1":
                # use mat[1, 1] for transverse scans, R0 factor is optional
                lorentz_factor = r0 * np.sqrt(det) / np.sqrt(mat[1, 1]) / np.sqrt(2 * np.pi)
            elif ax == "th2th":
                lorentz_factor = r0 * np.sqrt(det) / np.sqrt(mat[0, 0]) / np.sqrt(2 * np.pi)
            else:
                raise ValueError("axis not defined. Only supporting export s1 and the2th scans currently.")
            # ====================================================
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
            target = save_to_file if overwrite else VERITAS._next_version_path(save_to_file)
            with open(target, "w") as f:
                f.write(f"{title}\n")
                f.write("(3i5,2f8.2,i4,3f8.2)\n")
                f.write(f"{wavelength}  0   0\n")

                for (h, k, l), intensity, err in export:
                    f.write(f"{int(round(h)):5d}{int(round(k)):5d}{int(round(l)):5d}{intensity:8.2f}{err:8.2f}   1\n")
                f.close()
        return export
