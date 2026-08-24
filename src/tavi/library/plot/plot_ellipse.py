"""Handle plotting a 2d ellipse resolution matrix."""

from dataclasses import dataclass
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.library.experiment.utilities import sig2fwhm


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

        for entry in self.entries:
            entry.patch.set_transform(shear + ax.transData)
            ax.add_patch(entry.patch)
            x0, y0 = entry.origin
            xlo, ylo = x0 - pad * entry.x_extent, y0 - pad * entry.y_extent
            xhi, yhi = x0 + pad * entry.x_extent, y0 + pad * entry.y_extent

        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)

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
            ellipses: Sequence of (idx, peak, ellipse_co, axes_angle, coh_para,
                coh_perp, ellipse_incoh, incoh_para, incoh_perp) tuples. ``idx``
                and ``peak`` provide the scan number and hkl for the subplot
                title, ``ellipse_co`` / ``ellipse_incoh`` are the 2D resolution
                matrices drawn as a solid and a dotted ellipse respectively,
                ``axes_angle`` is the skew angle (degrees) between the plot
                axes, and ``coh_para`` / ``coh_perp`` / ``incoh_para`` /
                ``incoh_perp`` are the FWHMs shown in the title.
            xlabel: Custom x-axis label applied to each subplot. Left unlabeled if None.
            ylabel: Custom y-axis label applied to each subplot. Left unlabeled if None.

        """
        n = len(ellipses)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig = plt.figure(figsize=(4 * ncols, 4 * nrows))
        for i, (
            idx,
            peak,
            ellipse_co,
            axes_angle,
            coh_para,
            coh_perp,
            ellipse_incoh,
            incoh_para,
            incoh_perp,
        ) in enumerate(ellipses, start=1):
            ax = fig.add_subplot(nrows, ncols, i, axes_class=Axes, grid_helper=grid_helper(axes_angle))
            ax.grid(True)
            h, k, l = peak.hkl
            ax.set_title(
                f"{idx}, ({h:.2f} {k:.2f} {l:.2f})"
                f"\n coh_{xlabel} = {coh_para:.3f}, coh_{ylabel} = {coh_perp:.3f}"
                f"\n incoh_{xlabel} = {incoh_para:.3f}, incoh_{ylabel} = {incoh_perp:.3f}"
            )
            p = cls(axes_angle=axes_angle)
            p.add_ellipse(ellipse_co)
            p.add_ellipse(ellipse_incoh, linestyle=":")
            p.plot(ax=ax, show=False)
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.show()
