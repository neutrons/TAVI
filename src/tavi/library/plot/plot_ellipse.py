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


class Plot:
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
