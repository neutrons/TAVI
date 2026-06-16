"""General utilities for tas related functions and classes."""

from enum import Enum
from typing import Any, List, Optional, Tuple

import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.peak import DataPoint
from tavi.library.fit.fit import FitPackage, ModelName
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
        plugin: Any = None,
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
            plugin: Optional instrument-specific plugin used for data reduction.

        """
        self.instrument = instrument
        self.experiment = experiment
        self.sample = sample
        self.ub_convention = ub_convention
        self.ub_algorithm = UBAlgorithm(self.sample, self.instrument)
        self.plugin = plugin

    def ub(
        self,
        peaks: tuple[DataPoint, ...],
        scattering_plane: Optional[tuple[tuple, tuple]] = None,
        reset_ub: bool = False,
    ) -> np.ndarray:
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
        if reset_ub:
            self.sample.ol.UB = ub
        return ub

    def calculate_resolution(
        self, model: ResolutionModel = ResolutionModel.CooperNathans, axes: Optional[tuple] = None
    ) -> None:
        """Utilize Resolution class."""
        self.resolution = Resolution(
            model=model, instrument=self.instrument, sample=self.sample, experiment=self.experiment, axes=axes
        )

    def browse(
        self,
        scan_list: list[int],
        show_fits: bool = True,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: List[Tuple] = [],
        with_resolution_bar: bool = False,
        show_components: bool = False,
        def_x: str = None,
        def_y: str = None,
        xlim: float | list[float] = None,
        ylim: float | list[float] = None,
    ) -> None:
        """
        Browse scan with options to show resolution bar.

        If resolution bar is set, show_fits must be true as reso bar's positions are
        determined by fit center. When show_components is set, each fitted component
        (peak, linear background, ...) is also plotted separately. A scalar xlim/ylim
        pads the axis symmetrically around the data range (e.g. 1.1 widens the x-axis
        to 1.1x the data span about its midpoint); a [min, max] list sets it directly.
        """
        resolution_bar_4d = None
        if with_resolution_bar:
            resolution_bar_4d = self.resolution_bar(scan_list)
            # if with_resolution_bar is turned on, return fit_resutls and res_4d to be prepared for
            # intensity export.
            return browse_scans(
                self.experiment,
                scan_list,
                show_fits,
                fit_package,
                model_dict,
                resolution_bar_4d,
                show_components,
                def_x,
                def_y,
                xlim,
                ylim,
            )
        else:
            browse_scans(
                self.experiment,
                scan_list,
                show_fits,
                fit_package,
                model_dict,
                resolution_bar_4d,
                show_components,
                def_x,
                def_y,
                xlim,
                ylim,
            )

    def browse_resolution_ellipse(
        self,
        scan_list: list[int],
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: List[Tuple] = [(ModelName.Gaussian, dict(guess=True))],
    ) -> None:
        """
        Plot the resolution ellipse for each scan in a grid of subplots.

        Args:
            scan_list: Scan numbers to plot.
            xlabel: Custom x-axis label for each subplot. Left unlabeled if None.
            ylabel: Custom y-axis label for each subplot. Left unlabeled if None.
            fit_package: Fitting backend used to locate each peak center.
            model_dict: Models passed to the fit when locating peak centers.

        """
        if isinstance(self.experiment.loader, ORNLSpiceLoader):
            peaks = [
                dp
                for i in scan_list
                for dp in self.experiment.get_peak_center(
                    dict(scan_num=i), fit_package=fit_package, model_dict=model_dict
                )
            ]
        else:
            raise ValueError("Data format not implemented yet.")

        self.calculate_resolution()
        ellipses = []
        for idx, peak in zip(scan_list, peaks):
            res_4d, r0 = self.resolution.get_resolution(hkl=peak.hkl, ei=peak.ei, ef=peak.ef, rot_mat=None)
            ellipse, axes_angle = self.resolution.get_ellipse(res_mat=res_4d, ellipse_axes=(0, 1))
            coh_para = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=0)
            coh_perp = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=1)
            ellipses.append((idx, peak, ellipse, axes_angle, coh_para, coh_perp))

        PlotResolution.plot_resolution_ellipse(ellipses, xlabel=xlabel, ylabel=ylabel)

    # TODO: implement model_dict instead of just guessing, implement projection, axes parameters.
    def resolution_bar(self, scan_list: list[int], ax: int = 0) -> tuple[list[float], list]:
        """Compute the coherent FWHM resolution bar for each scan."""
        self.calculate_resolution()
        resolution_bar = []
        res_4ds = []
        if isinstance(self.experiment.loader, ORNLSpiceLoader):
            peaks = [
                dp
                for i in scan_list
                for dp in self.experiment.get_peak_center(
                    {"scan_num": i}, FitPackage.lmfit, [(ModelName.Gaussian, dict(guess=True))]
                )
            ]
        else:
            raise ValueError("Data format not implemented yet.")
        for peak in peaks:
            res_4d, r0 = self.resolution.get_resolution(hkl=peak.hkl, ei=peak.ei, ef=peak.ef, rot_mat=None)
            coh = ResolutionEllipsoid(res_4d, axes=None).coh_fwhm(axis=ax)
            resolution_bar.append(coh)
            res_4ds.append((res_4d, r0))
        return resolution_bar, res_4ds
