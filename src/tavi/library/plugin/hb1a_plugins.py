"""Hb1a specific plugins to reduce diffraction data."""

from pathlib import Path
from typing import Optional

import numpy as np


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
    def export_intensity(
        hkls: list,
        fit_results: list,
        res_4ds: list,
        ax: str,
        save_to_file: Optional[str],
        wavelength: float = 2.37815,
        overwrite: bool = True,
    ) -> list:
        """
        Export the intensity data to a .int file for refinement.

        Args:
            hkls: List of (h, k, l) per peak.
            fit_results: Fit result per peak, providing amplitude and amplitude_err.
            res_4ds: (resolution matrix, r0) per peak.
            ax: Scan axis; "s1" for transverse scans.
            save_to_file: Output file path. No file is written if None.
            wavelength: Neutron wavelength written to the file header.
            overwrite: If True (default), write to save_to_file, replacing it if it
                exists. If False, write to a new file with "_<n>" appended before the
                suffix, where n is the next unused version number on disk.

        """
        export: list = []
        for hkl, fit_result, res_4d in zip(hkls, fit_results, res_4ds):
            mat, r0 = res_4d[0], res_4d[1]
            # ====================================================
            # resolution calculated here
            det = np.abs(mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0])
            if ax == "s1":
                # use mat[1, 1] for transverse scans, R0 factor is optional
                lorentz_factor = r0 * np.sqrt(det) / np.sqrt(mat[1, 1]) / np.sqrt(2 * np.pi)
            # ====================================================
            peak = fit_result.peak
            intensity = peak.values["amplitude"] / lorentz_factor
            err = peak.errors["amplitude"] / lorentz_factor
            export.append((hkl, intensity, err))
        if save_to_file:
            target = save_to_file if overwrite else VERITAS._next_version_path(save_to_file)
            with open(target, "w") as f:
                f.write("Single crystal data of NdBi (hb1a)\n")
                f.write("(3i5,2f8.2,i4,3f8.2)\n")
                f.write(f"{wavelength}  0   0\n")

                for (h, k, l), intensity, err in export:
                    f.write(f"{int(round(h)):5d}{int(round(k)):5d}{int(round(l)):5d}{intensity:8.2f}{err:8.2f}   1\n")
                f.close()
        return export
