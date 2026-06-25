"""Handle plotting a 2d ellipse resolution matrix."""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.utilities import sig2fwhm
from tavi.library.fit.fit import FitPackage


def browse_scans(
    experiment: Experiment,
    scan_list: list[int],
    show_fits: bool = True,
    fit_package: FitPackage = FitPackage.lmfit,
    model_dict: List[Tuple] = [],
    resolution_bars: Optional[list[float]] = None,
    show_components: bool = False,
    def_x: str = None,
    def_y: str = None,
    xlim: Optional[float | list[float]] = None,
    ylim: Optional[float | list[float]] = None,
) -> None:
    """
    Plot a grid of scans, optionally with Gaussian fits and resolution bars.

    Args:
        experiment: Experiment holding the loaded scans.
        scan_list: Scan numbers to plot, one subplot each.
        show_fits: If True, overlay a Gaussian fit on each scan.
        fit_package: Fitting backend used to fit each scan.
        model_dict: Models passed to the fit for each scan.
        resolution_bars: Optional per-scan coherent FWHM (in q), aligned with
            scan_list. When given, the x-axis switches to del_q and a
            resolution bar is drawn on each subplot at the fitted peak's
            half-maximum. Requires show_fits=True.
        show_components: If True, also plot each fitted component (e.g. the peak
            and the linear background) separately as dashed lines, in addition to
            the composite fit. Requires show_fits=True.
        def_x: Default x-axis motor/column name. Falls back to the scan metadata
            default when not given.
        def_y: Default y-axis detector/column name. Falls back to the scan metadata
            default when not given.
        xlim: A scalar pads the x-axis symmetrically around the data range (e.g. 1.1
            widens it to 1.1x the data span about its midpoint); a [min, max] list
            sets the x range directly.
        ylim: A scalar pads the y-axis symmetrically around the data range (e.g. 1.2
            widens it to 1.2x the data span about its midpoint); a [min, max] list
            sets the y range directly.

    """
    from tavi.library.fit.fit import Fit

    show_resolution_bar = resolution_bars is not None
    if show_resolution_bar and not show_fits:
        raise ValueError("resolution_bars requires show_fits=True (the bar is placed using the fit).")

    n = len(scan_list)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    fit = Fit(package=fit_package)
    fit_results = []
    hkls = []
    bars, res_mat_4d = resolution_bars if show_resolution_bar else ([0] * n, [0] * n)
    for ax, num, coh in zip(axes_flat, scan_list, bars):
        scan = experiment.get_data_from_scan_number(dict(scan_num=num))
        if not def_x:
            def_x = scan.metadata.def_x
        if not def_y:
            def_y = scan.metadata.def_y

        if def_x[0].isdigit():
            def_x = "_" + def_x
        if def_y[0].isdigit():
            def_y = "_" + def_y

        # With a resolution bar, plot against del_q (the bar width is in q);
        # otherwise plot against the raw default-x motor.
        if show_resolution_bar:
            x = np.asarray(experiment.get_delta_q(dict(scan_num=num)))
            xlabel = f"del_q({def_x})"
        else:
            x = np.asarray(scan.data.data[def_x])
            xlabel = def_x
        y = np.asarray(scan.data.data[def_y])

        ax.errorbar(x, y, yerr=np.sqrt(y), fmt="o")

        if show_fits:
            fit_result = fit.fit(x, y, model_dict)
            fit_results.append(fit_result)
            peaks = fit_result.peaks
            x_fine = np.linspace(x.min(), x.max(), 300)
            fwhm_label = ", ".join(f"{peak.values['fwhm']:.4g}" for peak in peaks)
            ax.plot(x_fine, fit_result.raw.eval(x=x_fine), label=f"fit={fwhm_label}")

            # Optionally overlay each component (peak, linear background, ...) separately.
            if show_components and len(fit_result.components) > 1:
                for prefix, y_comp in fit_result.raw.eval_components(x=x_fine).items():
                    ax.plot(x_fine, y_comp, "--", lw=1, label=prefix.rstrip("_") or "component")

            if show_resolution_bar:
                # Resolution bar: horizontal line of width = coherent FWHM (coh) in q,
                # centered on each fitted peak and drawn at the peak's half-maximum.
                # The half-maximum sits at the background plus half the peak height, so
                # any linear (or other non-peak) background is added beneath the bar.
                def is_background(key: str) -> bool:
                    # eval_components keys a non-empty prefix by the prefix itself and an
                    # empty-prefix component by its function name (e.g. "linear"). A
                    # component is background when it has no fitted center.
                    component = fit_result.components.get(key) or fit_result.components.get("")
                    return component is not None and "center" not in component.values

                for idx, peak in enumerate(peaks):
                    center = peak.values["center"]
                    components = fit_result.raw.eval_components(x=center)
                    background = sum(float(value) for key, value in components.items() if is_background(key))
                    bar_y = background + peak.values["height"] / 2
                    ax.errorbar(
                        center,
                        bar_y,
                        xerr=coh / 2,
                        color="red",
                        capsize=4,
                        lw=2,
                        label=f"reso={coh:.4g}" if idx == 0 else None,
                    )
            ax.legend(fontsize=8, loc="upper right")
        hkls.append(experiment.get_hkl(dict(scan_num=num)))
        ax.set_title(f"{num}, {experiment.get_hkl(dict(scan_num=num))}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(def_y)

        # A scalar scales the axis symmetrically around the data midpoint (e.g. 1.1
        # widens to 1.1x the data span); a [min, max] sequence sets the range directly.
        if xlim is not None:
            if np.isscalar(xlim):
                x_mid, x_half = (x.max() + x.min()) / 2, (x.max() - x.min()) / 2
                ax.set_xlim(x_mid - xlim * x_half, x_mid + xlim * x_half)
            else:
                ax.set_xlim(xlim[0], xlim[1])
        if ylim is not None:
            if np.isscalar(ylim):
                y_mid, y_half = (y.max() + y.min()) / 2, (y.max() - y.min()) / 2
                ax.set_ylim(y_mid - ylim * y_half, y_mid + ylim * y_half)
            else:
                ax.set_ylim(ylim[0], ylim[1])

    # Hide any leftover empty axes in the grid.
    for ax in axes_flat[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    plt.show()
    if resolution_bars is not None:
        return hkls, fit_results, res_mat_4d


def grid_helper(angle: float, nbins: tuple[int, int] = (5, 5)) -> GridHelperCurveLinear:
    """Build a curve linear grid helper that skews axes by ``angle`` degrees."""
    # Forward: (x, y) -> (x + y/tan(angle), y)
    # Equivalent to a horizontal skew by (90 - angle) degrees in the Affine2D().skew_deg function,
    # which applies (x, y) -> (x + y*tran(angle), y) transformation.
    return GridHelperCurveLinear(
        Affine2D().skew_deg(90 - angle, 0),
        grid_locator1=MaxNLocator(nbins=nbins[0], steps=[1, 2, 5]),
        grid_locator2=MaxNLocator(nbins=nbins[1], steps=[1, 2, 5]),
    )


@dataclass
class EllipseEntry:
    """One ellipse to draw + the extent it occupies in data coords."""

    patch: Ellipse
    x_extent: float
    y_extent: float
    origin: tuple[float, float]


class PlotResolution:
    """Accumulate 2D resolution ellipses and render them on a skewed grid."""

    def __init__(self, axes_angle: float) -> None:
        """Initialize with the skew angle (degrees) between the two plot axes."""
        self.axes_angle = axes_angle
        self.entries: list[EllipseEntry] = []

    @staticmethod
    def create_ellipse(
        mat: np.ndarray,
        origin: tuple[float, float] = (0.0, 0.0),
        **kwargs: object,
    ) -> EllipseEntry:
        """Build an Ellipse patch (with FWHM dimensions) and its extent."""
        eigvals, eigvecs = np.linalg.eigh(mat)
        fwhm = sig2fwhm / np.sqrt(np.abs(eigvals))  # sig2fwhm is already "full" width, no need to * 2.
        angle_rad = np.arctan2(eigvecs[1, 0], eigvecs[0, 0])
        angle_deg = np.degrees(angle_rad)

        patch_kwargs = {"fill": False, "edgecolor": "0", "lw": 1.5, **kwargs}
        patch = Ellipse(
            xy=origin,
            width=fwhm[0],
            height=fwhm[1],
            angle=angle_deg,
            **patch_kwargs,
        )

        x_extent = np.sqrt((fwhm[0] / 2 * np.cos(angle_rad)) ** 2 + (fwhm[1] / 2 * np.sin(angle_rad)) ** 2)
        y_extent = np.sqrt((fwhm[0] / 2 * np.sin(angle_rad)) ** 2 + (fwhm[1] / 2 * np.cos(angle_rad)) ** 2)
        return EllipseEntry(patch=patch, x_extent=x_extent, y_extent=y_extent, origin=origin)

    def add_ellipse(
        self,
        mat: np.ndarray,
        origin: tuple[float, float] = (0.0, 0.0),
        **kwargs: object,
    ) -> EllipseEntry:
        """Create an ellipse and append it to the draw queue. Returns the entry."""
        entry = self.create_ellipse(mat, origin=origin, **kwargs)
        self.entries.append(entry)
        return entry

    def plot(self, ax: Optional[Axes] = None, pad: float = 1.1, show: bool = True) -> Axes:
        """Draw every added ellipse on a skewed-grid axes."""
        if not self.entries:
            raise RuntimeError("No ellipses added. Call add_ellipse(...) first.")

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(
                1,
                1,
                1,
                axes_class=Axes,
                grid_helper=grid_helper(self.axes_angle),
            )
            ax.grid(True)

        shear = Affine2D().skew_deg(90 - self.axes_angle, 0)

        x_min, x_max = np.inf, -np.inf
        y_min, y_max = np.inf, -np.inf

        for entry in self.entries:
            entry.patch.set_transform(shear + ax.transData)
            ax.add_patch(entry.patch)
            x0, y0 = entry.origin
            x_min = min(x_min, x0 - pad * entry.x_extent)
            x_max = max(x_max, x0 + pad * entry.x_extent)
            y_min = min(y_min, y0 - pad * entry.y_extent)
            y_max = max(y_max, y0 + pad * entry.y_extent)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        # Show legend if any patch has a non-private label
        if any(e.patch.get_label() and not e.patch.get_label().startswith("_") for e in self.entries):
            ax.legend()

        plt.tight_layout()
        if show:
            plt.show()
        return ax

    @classmethod
    def plot_resolution_ellipse(
        cls,
        ellipses: list[tuple],
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
    ) -> None:
        """
        Lay out a grid of resolution ellipses, one subplot per peak.

        Args:
            ellipses: Sequence of (peak, mat, axes_angle, coh_para, coh_perp)
                tuples. ``peak`` provides the hkl for the subplot title, ``mat``
                is the 2D resolution matrix to draw, ``axes_angle`` is the skew
                angle (degrees) between the plot axes, and ``coh_para`` /
                ``coh_perp`` are the coherent FWHMs shown in the title.
            xlabel: Custom x-axis label applied to each subplot. Left unlabeled if None.
            ylabel: Custom y-axis label applied to each subplot. Left unlabeled if None.

        """
        n = len(ellipses)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig = plt.figure(figsize=(4 * ncols, 4 * nrows))
        for i, (idx, peak, mat, angle, coh_para, coh_perp) in enumerate(ellipses, start=1):
            ax = fig.add_subplot(nrows, ncols, i, axes_class=Axes, grid_helper=grid_helper(angle))
            ax.grid(True)
            h, k, l = peak.hkl
            ax.set_title(f"{idx}, ({h:g} {k:g} {l:g}), fwhm_s1 = {coh_perp:.4f}, \n fwhm_th2th = {coh_para:.4f}")
            p = cls(axes_angle=angle)
            p.add_ellipse(mat)
            p.plot(ax=ax, show=False)
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.show()
