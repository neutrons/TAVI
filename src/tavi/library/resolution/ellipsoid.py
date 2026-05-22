"""Manage 4D ellipsoid and its projection."""

from typing import Optional, Tuple

import numpy as np

from tavi.library.experiment.utilities import quadric_proj, sig2fwhm


class ResoEllipsoid:
    """4D ellipsoid and its projection."""

    def __init__(
        self,
        res_mat: Optional[np.ndarray] = None,
        r0: Optional[float] = None,
        axes: Optional[Tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"),
    ) -> None:
        """Init."""
        self.res_mat = res_mat
        self.r0 = r0
        self.axes = axes
        self.angles = (90, 90, 90)

    def project_to_frame(self, r_mat: np.ndarray, psi: float, ub: np.ndarray) -> np.ndarray:
        """
        Project the res_mat to a certain frame in hkle.

        Args:
            r_mat: the rotation matrix in UB formalism.
            psi: angle between Q and ki(z).
            ub: UB matrix of the sample's oriented lattice.

        """
        mat_lab_to_local = np.array(
            [
                [np.sin(psi), 0, np.cos(psi)],
                [np.cos(psi), 0, -np.sin(psi)],
                [0, 1, 0],
            ]
        )
        # hkle frame
        if self.axes == ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"):
            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = 2 * np.pi * mat_lab_to_local @ r_mat @ ub
            res_mat_proj = conv_mat_4d.T @ self.res_mat @ conv_mat_4d

        else:  # if project to other frames, multiply further by matrix W
            w1, w2, w3 = np.array([item for item in self.axes if isinstance(item, tuple) and len(item) == 3])
            # Eq.71
            mat_w = np.array([w1, w2, w3]).T
            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = 2 * np.pi * mat_lab_to_local @ r_mat @ ub @ mat_w
            res_mat_proj = conv_mat_4d.T @ self.res_mat @ conv_mat_4d

        return res_mat_proj, self.r0

    def coh_fwhm(self, res_mat: np.ndarray, axis: int = None) -> float:
        """
        Coherent FWHM.

        Make a cut.
        """
        idx = int(axis)
        return sig2fwhm / np.sqrt(res_mat[idx, idx])

    def incoh_fwhm(self, res_mat: np.ndarray, axis: int = None) -> float:
        """
        Incoherent FWHM.

        Integrate all 3 axes.
        """
        idx = int(axis)
        res_mat = res_mat
        for i in (3, 2, 1, 0):
            if not i == idx:
                res_mat = quadric_proj(res_mat, i)
        return sig2fwhm / np.sqrt(np.abs(res_mat[0, 0]))
