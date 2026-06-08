"""General utilities for tas related functions and classes."""

from enum import Enum
from typing import Optional

import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.peak import DataPoint
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument
from tavi.library.plot.plot_ellipse import PlotResolution, browse_scans
from tavi.library.resolution.ellipsoid import ResolutionEllipsoid
from tavi.library.resolution.resolution import Resolution, ResolutionModel
from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader
from tavi.library.ubalgorithm.ub import UBAlgorithm


class UBConvention(Enum):
    """Convention used to define the UB matrix."""

    Spice = "Spice convention."
    Mantid = "Mantid convention."


class TAS:
    """
    Triple-axis class. Main function for tavi library.

    Next step is implement resolution calculation from tavi.library.resolution module.
    """

    def __init__(
        self,
        instrument: Instrument,
        sample: Sample,
        experiment: Optional[Experiment] = None,
        ub_convention: UBConvention = UBConvention.Spice,
    ) -> None:
        """
        Initialize triple axis.

        Args:
            instrument: instrument name. Currently placeholder but will be used to load instrument
            configuration later.
            goni: Goniometer.
            sample: Sample, need to construct OrientedLattice first. See sample.py.
            resolution: Resolution method to use for resolution calculations.
            experiment: Experiment. Handles extracting peak center, ei, ef etc from exp data.
            ub_convention: Convention used to define the UB matrix. See UBConvention.

        """
        self.instrument = instrument
        self.experiment = experiment
        self.sample = sample
        self.ub_convention = ub_convention
        self.ub_algorithm = UBAlgorithm(self.sample, self.instrument)

    def ub(self, peaks: tuple[DataPoint, ...], scattering_plane: Optional[tuple[tuple, tuple]] = None) -> np.ndarray:
        """Calculate ub matrix from 1,2,3 or multiple peaks."""
        if len(peaks) == 1 and scattering_plane is None:
            raise ValueError("Scattering plane cannot be None with only 1 peak specified.")

        match len(peaks):
            case 1:
                # call find u from one peak and scattering plane
                peak = peaks[0]
                u = self.ub_algorithm.find_u_from_one_peak_and_scattering_plane(
                    peak=peak, scattering_plane=scattering_plane, ei=peak.ei, ef=peak.ef
                )
                ub = u @ self.sample.ol.B
            case 2:
                # call find u from 2 peaks
                u = self.ub_algorithm.find_u_from_two_peaks(peaks)
                ub = u @ self.sample.ol.B
            case 3:
                # call find ub from 3 peaks
                ub = self.ub_algorithm.find_ub_from_three_peaks(peaks)
            case _:
                ub = self.ub_algorithm.find_ub_from_multiple_peaks(peaks)
        self.sample.ol.UB = ub
        return ub

    def calculate_resolution(
        self, model: ResolutionModel = ResolutionModel.CooperNathans, axes: Optional[tuple] = None
    ) -> None:
        """Utilize Resolution class."""
        self.resolution = Resolution(
            model=model, instrument=self.instrument, sample=self.sample, experiment=self.experiment, axes=axes
        )

    def browse(self, scan_list: list[int], show_fits=True, with_reso_bar=False):
        """
        Browse scan with options to show resolution bar. If resolution bar is set,
        show_fits must be true as reso bar's poitions are deteremined by fit center.
        """
        reso_bar = None
        if with_reso_bar:
            reso_bar = self.reso_bar(scan_list)
        browse_scans(self.experiment, scan_list, show_fits, reso_bar)

    def browse_resolution_ellipse(self, scan_list: list[int]) -> None:
        """Plot the resolution ellipse for each scan in a grid of subplots."""
        if isinstance(self.experiment.loader, ORNLSpiceLoader):
            peaks = [self.experiment.create_peaks(dict(scan_num=i)) for i in scan_list]
        else:
            raise ValueError("Data format not implemented yet.")

        self.calculate_resolution()
        ellipses = []
        for peak in peaks:
            res_4d, r0 = self.resolution.get_resolution(hkl=peak.hkl, ei=peak.ei, ef=peak.ef, rot_mat=None)
            ellipse, axes_angle = self.resolution.get_ellipse(res_mat=res_4d, ellipse_axes=(0, 1))
            coh_para = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=0)
            coh_perp = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=1)
            ellipses.append((peak, ellipse, axes_angle, coh_para, coh_perp))

        PlotResolution.plot_reso_ellipse(ellipses)

    def reso_bar(self, scan_list: list[int]):
        self.calculate_resolution()
        reso_bar = []
        if isinstance(self.experiment.loader, ORNLSpiceLoader):
            peaks = [self.experiment.create_peaks(dict(scan_num=i)) for i in scan_list]
        else:
            raise ValueError("Data format not implemented yet.")
        for peak in peaks:
            res_4d, r0 = self.resolution.get_resolution(hkl=peak.hkl, ei=peak.ei, ef=peak.ef, rot_mat=None)
            coh = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=0)
            reso_bar.append(coh)
        return reso_bar
