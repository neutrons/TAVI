"""Cooper-Nathans Method."""

from typing import Tuple

import numpy as np

from tavi.library.component.analyzer import Analyzer
from tavi.library.component.collimators import Collimators
from tavi.library.component.monochromater import Monochromator
from tavi.library.experiment.utilities import SE2K, ksq2eng, rotation_matrix_2d, sig2fwhm
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument


class CooperNathans:
    """
    Cooper-Nathans method.

    Follows Popvici 1975 paper. Calculate matrix G, F, C, A, B.
    """

    NUM_MONO = 1
    NUM_ANA = 1
    IDX_MONO = {"0H": 0, "0V": 1}
    IDX_ANA = {"0H": 2, "0V": 3}

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLLS = {"0H": 0, "0V": 2, "1H": 1, "1V": 3, "2H": 4, "2V": 6, "3H": 5, "3V": 7}

    def mat_g(self, collimators: Collimators) -> np.ndarray:
        """Matrix G, 8 x 8 diagonal matrix about the horizontal and vertical divergence of 4 collimators."""
        mat_g = np.zeros((self.NUM_COLLS * 2, self.NUM_COLLS * 2))

        mat_g[self.IDX_COLLS["0H"], self.IDX_COLLS["0H"]] = 1.0 / collimators.pre_mono_h**2
        mat_g[self.IDX_COLLS["0V"], self.IDX_COLLS["0V"]] = 1.0 / collimators.pre_mono_v**2

        mat_g[self.IDX_COLLS["1H"], self.IDX_COLLS["1H"]] = 1.0 / collimators.pre_sample_h**2
        mat_g[self.IDX_COLLS["1V"], self.IDX_COLLS["1V"]] = 1.0 / collimators.pre_sample_v**2

        mat_g[self.IDX_COLLS["2H"], self.IDX_COLLS["2H"]] = 1.0 / collimators.post_sample_h**2
        mat_g[self.IDX_COLLS["2V"], self.IDX_COLLS["2V"]] = 1.0 / collimators.post_sample_v**2

        mat_g[self.IDX_COLLS["3H"], self.IDX_COLLS["3H"]] = 1.0 / collimators.post_ana_h**2
        mat_g[self.IDX_COLLS["3V"], self.IDX_COLLS["3V"]] = 1.0 / collimators.post_ana_v**2

        return mat_g

    def mat_f(self, monochromator: Monochromator, analyzer: Analyzer) -> np.ndarray:
        """
        Matrix F.

        A 4 x 4 diagonal matrix considering divergence of monochromator and analyzer.
        """
        mat_f = np.zeros(((self.NUM_MONO + self.NUM_ANA) * 2, (self.NUM_MONO + self.NUM_ANA) * 2))
        mat_f[self.IDX_MONO["0H"], self.IDX_MONO["0H"]] = 1.0 / monochromator.mosaic_h**2
        mat_f[self.IDX_MONO["0V"], self.IDX_MONO["0V"]] = 1.0 / monochromator.mosaic_v**2
        mat_f[self.IDX_ANA["0H"], self.IDX_ANA["0H"]] = 1.0 / analyzer.mosaic_h**2
        mat_f[self.IDX_ANA["0V"], self.IDX_ANA["0V"]] = 1.0 / analyzer.mosaic_v**2
        return mat_f

    def mat_a(self, ki: float, kf: float, theta_m: float, theta_a: float) -> np.ndarray:
        """
        Calculate Matrix A. 6 x 8 matrix, Y = AU. Transform collmiator's divergence to ki-kf.

        Args:
            ki: in meV.
            kf: in meV.
            theta_m: in radians.
            theta_a: in radians.

        """
        mat_a = np.zeros((6, 2 * self.NUM_COLLS))
        mat_a[0, self.IDX_COLLS["0H"]] = 0.5 * ki / np.tan(theta_m)
        mat_a[0, self.IDX_COLLS["1H"]] = -0.5 * ki / np.tan(theta_m)
        mat_a[1, self.IDX_COLLS["1H"]] = ki
        mat_a[2, self.IDX_COLLS["1V"]] = ki

        mat_a[3, self.IDX_COLLS["2H"]] = 0.5 * kf / np.tan(theta_a)
        mat_a[3, self.IDX_COLLS["3H"]] = -0.5 * kf / np.tan(theta_a)
        mat_a[4, self.IDX_COLLS["2H"]] = kf
        mat_a[5, self.IDX_COLLS["2V"]] = kf
        return mat_a

    def mat_b(self, ki: float, kf: float, phi: float, two_theta: float) -> np.ndarray:
        """
        Calculate Matrix B. 4 x 6 matrix, X = BY. Transform from ki-kf to q-frame.

        Args:
            ki: in meV.
            kf: in meV.
            phi: in radians.
            two_theta: in radians.

        """
        mat_b = np.zeros((4, 6))
        mat_b[0:3, 0:3] = rotation_matrix_2d(phi)
        mat_b[0:3, 3:6] = rotation_matrix_2d(phi - two_theta) * (-1)
        mat_b[3, 0] = 2 * ksq2eng * ki
        mat_b[3, 3] = -2 * ksq2eng * kf
        return mat_b

    def mat_c(self, theta_m: float, theta_a: float) -> np.ndarray:
        """
        Matrix C. 4 x 8 matrix. Constrinat between mono/ana mosaic and collimator divergence.

        Args:
            theta_m: in radians.
            theta_a: in radians.

        """
        mat_c = np.zeros(((self.NUM_MONO + self.NUM_ANA) * 2, self.NUM_COLLS * 2))
        mat_c[self.IDX_MONO["0H"], self.IDX_COLLS["0H"]] = 0.5
        mat_c[self.IDX_MONO["0H"], self.IDX_COLLS["1H"]] = 0.5
        mat_c[self.IDX_MONO["0V"], self.IDX_COLLS["0V"]] = 0.5 / np.sin(theta_m)
        mat_c[self.IDX_MONO["0V"], self.IDX_COLLS["1V"]] = -0.5 / np.sin(theta_m)

        mat_c[self.IDX_ANA["0H"], self.IDX_COLLS["2H"]] = 0.5
        mat_c[self.IDX_ANA["0H"], self.IDX_COLLS["3H"]] = 0.5
        mat_c[self.IDX_ANA["0V"], self.IDX_COLLS["2V"]] = 0.5 / np.sin(theta_a)
        mat_c[self.IDX_ANA["0V"], self.IDX_COLLS["3V"]] = -0.5 / np.sin(theta_a)
        return mat_c

    def resolution_matrix(
        self, instrument: Instrument, sample: Sample, hkl: Tuple[float, float, float], ei: float, ef: float
    ) -> tuple[np.ndarray, float]:
        """Calculate resolution matrix and normalization factor based on Popvic 1975."""
        # To be implemented: instrument will contain information about collimators, monochromator, analyzer, angles etc.

        q_norm = sample.ol.q_norm_from_hkl(hkl)
        ki, kf = SE2K(ei), SE2K(ef)
        # motor_angles = instrument.calculate_motor_angles(hkl=hkl, en=ei - ef)
        # psi = <ki to q>, always has the oppositie sign of s2
        psi = instrument.get_psi(hkl, sample, ei, ef)
        two_theta = instrument.get_two_theta(hkl, sample, ei, ef)
        theta_m = instrument.monochromater.theta_m(ei) * (
            1 if instrument.monochromater.sense == "+" else -1
        )  # set sense
        theta_a = instrument.analyzer.theta_a(ef) * (1 if instrument.analyzer.sense == "+" else -1)  # set sense

        mat_a = self.mat_a(ki, kf, theta_m, theta_a)
        mat_b = self.mat_b(ki, kf, psi, two_theta)
        mat_c = self.mat_c(theta_m, theta_a)
        mat_f = self.mat_f(instrument.monochromater, instrument.analyzer)
        mat_g = self.mat_g(instrument.collimators)

        mat_h = mat_c.T @ mat_f @ mat_c + mat_g
        mat_h_inv = np.linalg.inv(mat_h)
        mat_ba = mat_b @ mat_a
        mat_cov = mat_ba @ mat_h_inv @ mat_ba.T

        mat_cov[1, 1] += q_norm**2 * sample.mosaic_h**2
        mat_cov[2, 2] += q_norm**2 * sample.mosaic_v**2
        # times sig2fwhm^2 to convert from sigma to FWHM
        mat_reso = np.linalg.inv(mat_cov) * sig2fwhm**2

        # TODO check normalization factor
        # -------------------------------------------------------------------------
        # - if the instruments works in kf=const mode and the scans are counted for
        #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
        # - if the instrument works in ki=const mode the kf^3 factor is needed.

        # monochromator and analyzer reflectivity taken to be 1
        rm_theta = 1
        ra_theta = 1
        rm_k = rm_theta * ki**3 / np.tan(theta_m)
        ra_k = ra_theta * kf**3 / np.tan(theta_a)

        r0 = rm_k * ra_k
        r0 *= np.pi**2 / 4 / np.sin(theta_m) / np.sin(theta_a)
        r0 *= np.sqrt(np.linalg.det(mat_f) / np.linalg.det(mat_h))

        return (mat_reso, r0)
