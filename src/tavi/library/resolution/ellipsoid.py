"""Manage 4D ellipsoid and its projection."""

from typing import Any, Optional, Tuple

import numpy as np

from tavi.library.experiment.utilities import sig2fwhm
from tavi.library.geometry.sample import Sample


class ResoEllipsoid:
    """4D ellipsoid and its projection."""

    def __init__(
        self,
        instrument: Any,
        sample: Sample,
        hkle: Tuple[float, float, float, float],
        axes: Optional[Tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"),
        res_mat: Optional[np.ndarray] = None,
        r0: Optional[float] = None,
    ) -> None:
        """Init."""
        self.instrument = instrument
        *self.hkl, self.e = hkle
        self.axes = axes
        self.res_mat = res_mat
        self.r0 = r0

        self.q = sample.ol.q_norm_from_hkl(self.hkl)

    def project_to_frame(self, r_mat: np.ndarray) -> np.ndarray:
        """
        Project the res_mat to a certain frame in hkle.

        Args:
            r_mat: the rotation matrix in UB formalism.

        """
        ub = self.sample.ol.UB
        # angle between Q and ki(z)
        psi = self.instrument.get_psi(self.hkl, self.e)

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
            conv_mat_4d[0:3, 0:3] = 2 * np.pi @ mat_lab_to_local @ r_mat @ ub
            res_mat_proj = conv_mat_4d.T @ self.res_mat @ conv_mat_4d

        else:  # if project to other frames, multiply further by matrix W
            w1, w2, w3 = np.array([item for item in self.axes if isinstance(item, tuple) and len(item) == 3])
            # Eq.71
            mat_w = np.array([w1, w2, w3]).T
            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = 2 * np.pi @ mat_lab_to_local @ r_mat @ ub @ mat_w
            res_mat_proj = conv_mat_4d.T @ self.res_mat @ conv_mat_4d

            # update axes
            e_axis = self.axes.index("e")
            new_idx = [0, 1, 2]
            new_idx.insert(e_axis, 3)
            res_mat_proj = res_mat_proj[new_idx, :][:, new_idx]
        return res_mat_proj

    def quadric_proj(self, quadric: np.ndarray, idx: int) -> np.ndarray:
        """
        Project out a specific dimension.

        Implementation is equivalent of R.A.Robinson et al. Analytical Techniques for Instrument
        Design - Matrix Methods Eq. 13 (demonimeter should not be squared). Or Eck[14] 57. Calculation follows Takin's ellipse.h example.

        We can verify with:
            k = 0
            test = np.zeros((4,4))
            for i in range(4):
                for j in range(4):
                    if i == k or j == k:
                        continue
                    test[i][j] = quadric[i][j] - quadric[i][k]*quadric[j][k]/quadric[k][k].

        """
        atol = 1e-8

        if np.abs(quadric[idx, idx]) < atol:
            return np.delete(np.delete(quadric, idx, axis=0), idx, axis=1)

        # row/column along which to perform the orthogonal projection
        vec = 0.5 * (quadric[idx, :] + quadric[:, idx])  # symmetrise if not symmetric
        vec /= np.sqrt(quadric[idx, idx])  # normalise to indexed component
        proj_op = np.outer(vec, vec)  # projection operator
        ortho_proj = quadric - proj_op  # projected quadric
        return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)

    def get_ellipse(
        self, r_mat: np.ndarray, ellipse_axes: tuple[int, int] = (0, 1), PROJECTION: bool = False, ORIGIN: bool = True
    ) -> np.ndarray:
        """Generate 2D resolution ellipse."""
        res_2d = r_mat
        for idx in (3, 2, 1, 0):
            if idx not in ellipse_axes:
                # make an integral
                if PROJECTION:
                    res_2d = self.quadric_proj(res_2d, idx)
                # make a slice/cut
                else:
                    res_2d = np.delete(np.delete(res_2d, idx, axis=0), idx, axis=1)
        return res_2d

    def coh_fwhm(self, r_mat: np.ndarray, axis: int = None) -> float:
        """
        Coherent FWHM.

        Make a cut.
        """
        idx = int(axis)
        return sig2fwhm / np.sqrt(r_mat[idx, idx])

    def incoh_fwhm(self, r_mat: np.ndarray, axis: int = None) -> float:
        """
        Incoherent FWHM.

        Integrate all 3 axes.
        """
        idx = int(axis)
        res_mat = r_mat
        for i in (3, 2, 1, 0):
            if not i == idx:
                res_mat = self.quadric_proj(res_mat, i)
        return sig2fwhm / np.sqrt(np.abs(res_mat[0, 0]))
