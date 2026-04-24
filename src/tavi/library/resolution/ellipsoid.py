"""Manage 4D ellipsoid and its projection."""

import numpy as np
from typing import Any, Optional, Tuple

from old_tavi.old_tests.test_lattice_ub_algorithm import sample_info
from tavi.library.geometry.sample import Sample

class ResoEllipsoid:
    """4D ellipsoid and its projection."""

    def __init__(self, 
                 instrument: Any, 
                 sample:Sample,
                 hkle: Tuple[float, float, float, float], 
                 axes: Optional[Tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"),
                 res_mat: Optional[np.ndarray] = None,
                 r0: Optional[float] = None):
        self.instrument = instrument
        *self.hkl, self.e = hkle
        self.axes = axes
        self.res_mat = res_mat
        self.r0 = r0
    
    def project_to_frame(self, r_mat:np.ndarray) -> np.ndarray:
        """
        Project the res_mat to a certain frame.
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

        if self.axes == ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"):
            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = 2 * np.pi @ mat_lab_to_local @ r_mat @ ub
            res_mat_proj = conv_mat_4d.T @ self.res_mat @ conv_mat_4d
        else: # if project to other frames, multiply further by matrix W
            pass
        return res_mat_proj