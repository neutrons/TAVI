"""General utilities for tas related functions and classes."""

from enum import Enum
from typing import Any, List, Optional, Tuple

import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.peak import DataPoint
from tavi.library.fit import FitPackage, ModelName
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument
from tavi.library.plot.browser import browse_scans
from tavi.library.plot.plot_ellipse import PlotResolution
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

    def browse(
        self,
        scan_list: list[int],
        show_fits: bool = True,
        fit_package: FitPackage = FitPackage.lmfit,
        model_dict: List[Tuple] = [],
        show_resolution_bar: bool = False,
        show_components: bool = False,
        def_x: str = None,
        def_y: str = None,
        xlim: float | list[float] = None,
        ylim: float | list[float] = None,
        resolution_frame: str | tuple = "local",
        projection_axis: int = 0,
    ) -> None:
        """
        Browse scan with options to show resolution bar.

        If resolution bar is set, show_fits must be true as reso bar's positions are
        determined by fit center. When show_components is set, each fitted component
        (peak, linear background, ...) is also plotted separately. A scalar xlim/ylim
        pads the axis symmetrically around the data range (e.g. 1.1 widens the x-axis
        to 1.1x the data span about its midpoint); a [min, max] list sets it directly.

        Args:
            scan_list: Scan numbers to browse.
            show_fits: Overlay the fitted curve on each scan. Must be True when
                show_resolution_bar is set.
            fit_package: Fitting backend used to fit each scan.
            model_dict: Models passed to the fit (peak, background, ...).
            show_resolution_bar: Overlay a resolution bar whose position is taken
                from the fit center; returns fit results and the 4D resolution for
                intensity export.
            show_components: Plot each fitted component separately in addition to
                the total fit.
            def_x: Name of the motor/field to use for the x-axis. Defaults to the
                scan's default if None.
            def_y: Name of the detector/field to use for the y-axis. Defaults to the
                scan's default if None.
            xlim: x-axis limits. A scalar pads symmetrically about the data midpoint;
                a [min, max] list sets the range directly.
            ylim: y-axis limits, interpreted like xlim.
            resolution_frame: Frame in which the resolution bar is computed.
                ``"local"`` for the local Q frame, or a 4-tuple of projection
                vectors plus ``"e"`` for an hkle frame, e.g.
                ``((1, 0, 0), (0, 1, 0), (0, 0, 1), "e")``.
            projection_axis: Axis of resolution_frame the FWHM is taken along.
                0/1/2 for the momentum axes, 3 for energy.

        """
        resolution_bar_4d = None
        if show_resolution_bar:
            # if show_resolution_bar is turned on, return fit_resutls and res_4d to be prepared for
            # intensity export.
            resolution_bar_4d = self.resolution_bar(
                scan_list, resolution_frame=resolution_frame, projection_axis=projection_axis, model_dict=model_dict
            )
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

    def resolution_bar(
        self,
        scan_list: list[int],
        resolution_model: ResolutionModel = ResolutionModel.CooperNathans,
        resolution_frame: str | Tuple = "local",
        projection_axis: int = 0,
        model_dict: List[Tuple] = [(ModelName.Gaussian, dict(guess=True))],
    ) -> tuple[list[float], list]:
        """
        Compute the coherent FWHM resolution bar for each scan.

        Args:
            scan_list: Scan numbers to compute the resolution for.
            resolution_model: Resolution model, default Cooper-Nathans.
            resolution_frame: ``"local"`` for the local Q frame, or a 4-tuple of
                projection vectors plus ``"e"`` for an hkle frame, e.g.
                ``((1, 0, 0), (0, 1, 0), (0, 0, 1), "e")``.
            projection_axis: Axis of resolution_frame the FWHM is taken along.
                0/1/2 for the momentum axes, 3 for energy.
            model_dict: Models used to fit each scan when locating peak centers.

        Returns:
            ``(resolution_bar, res_4ds)``. Both are lists with one entry per
            scan, each entry a tuple with one item per fitted peak in that scan:
            the coherent FWHM, and the ``(res_mat, r0)`` pair respectively.

        """
        # This step just get the center data point, has nothing to do with calculations.
        if isinstance(self.experiment.loader, ORNLSpiceLoader):
            # Group the data points from one scan i into a tuple, so centers is aligned
            # one-entry-per-scan rather than flattened across all scans.
            centers = [
                tuple(self.experiment.get_closest_to_center_data_point({"scan_num": i}, FitPackage.lmfit, model_dict))
                for i in scan_list
            ]
        else:
            raise ValueError("Data format not implemented yet.")
        resolution_bar = []
        incoh_bar = []  # leave it here to allow plotting incoh resolution bar in the future
        res_4ds = []

        # check what kind of projection_axes we are using. We supply two options:
        # 1. if it's 0 or 1, this implies user is requesting resolution bar in local Q frame.
        # This projection_axes goes through ellipsoid.py->ResolutionEllipsoid->coh/incoh_fwhm.
        # 2. if it's a tuple ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e") or ((1, 1, 0), (1, -1, 0), (0, 0, 1), "e"),
        # this implies user is requesting resolution in hkle frame.
        self.resolution = Resolution(
            model=resolution_model,
            instrument=self.instrument,
            sample=self.sample,
            experiment=self.experiment,
            resolution_frame=resolution_frame,
        )
        if resolution_frame == "local":
            for center_group in centers:
                cohs, incohs, res_group = [], [], []
                for center in center_group:
                    res_4d, r0 = self.resolution.get_resolution(
                        hkl=center.hkl, ei=center.ei, ef=center.ef, rot_mat=None
                    )
                    resolution_ellipsoid = ResolutionEllipsoid(res_4d, resolution_frame=resolution_frame)
                    cohs.append(resolution_ellipsoid.coh_fwhm(projection_axis))
                    incohs.append(resolution_ellipsoid.incoh_fwhm(projection_axis))
                    res_group.append((res_4d, r0))
                resolution_bar.append(tuple(cohs))
                incoh_bar.append(tuple(incohs))
                res_4ds.append(tuple(res_group))
            print("incoherent FWHM = ", incoh_bar)
        elif isinstance(resolution_frame, tuple):
            for center_group in centers:
                cohs, incohs, res_group = [], [], []
                for center in center_group:
                    # We need to calculate rotation matrix now. Two ways to do this, from goniometer.r_mat
                    # or from oriented_lattice.rot_matrix_with_minimal_tilt.
                    try:
                        rot_mat = self.instrument.goni.r_mat(center.angles)
                    except Exception as err:
                        rot_mat = None
                        raise ValueError("Rotation matrix can not be calculated") from err
                    res_4d, r0 = self.resolution.get_resolution(
                        hkl=center.hkl, ei=center.ei, ef=center.ef, rot_mat=rot_mat
                    )
                    resolution_ellipsoid = ResolutionEllipsoid(res_4d, resolution_frame=resolution_frame)
                    cohs.append(resolution_ellipsoid.coh_fwhm(projection_axis))
                    incohs.append(resolution_ellipsoid.incoh_fwhm(projection_axis))
                    res_group.append((res_4d, r0))
                resolution_bar.append(tuple(cohs))
                incoh_bar.append(tuple(incohs))
                res_4ds.append(tuple(res_group))
            # printing for now. Later it will be massaged into visualization.
            print("incoherent FWHM = ", incoh_bar)
        else:
            raise ValueError("Resolution frame is not defined properly.")
        return resolution_bar, res_4ds  # , incoh_bar

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
                (i, dp)
                for i in scan_list
                for dp in self.experiment.get_closest_to_center_data_point(
                    dict(scan_num=i), fit_package=fit_package, model_dict=model_dict
                )
            ]
        else:
            raise ValueError("Data format not implemented yet.")

        self.resolution = Resolution(
            model=ResolutionModel.CooperNathans,
            instrument=self.instrument,
            sample=self.sample,
            experiment=self.experiment,
            resolution_frame="local",
        )
        ellipses = []
        for idx, peak in peaks:
            res_4d, r0 = self.resolution.get_resolution(hkl=peak.hkl, ei=peak.ei, ef=peak.ef, rot_mat=None)
            ellipse, axes_angle = self.resolution.get_ellipse(res_mat=res_4d, ellipse_axes=(0, 1))
            coh_para = ResolutionEllipsoid(res_4d, resolution_frame="local").coh_fwhm(projection_axis=0)
            coh_perp = ResolutionEllipsoid(res_4d, resolution_frame="local").coh_fwhm(projection_axis=1)
            ellipses.append((idx, peak, ellipse, axes_angle, coh_para, coh_perp))

        PlotResolution.plot_resolution_ellipse(ellipses, xlabel=xlabel, ylabel=ylabel)
