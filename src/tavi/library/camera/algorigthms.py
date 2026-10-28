"""Reduce Anger camera frames to hkl, filter them, and locate Bragg peaks."""

from math import ceil
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from scipy import ndimage

from tavi.library.experiment.utilities import rot_y, rot_z


class AngerCameraAlgo:
    """Convert Anger camera points to hkl and find the peaks in them."""

    def __init__(self, anger_camera_data: dict, wavelength: float) -> None:
        """Hold a parsed loader store and the incident wavelength in angstroms."""
        self.data = anger_camera_data
        self.wavelength = wavelength

    def get_scans(self, name: str) -> dict:
        """Return {filename: data} for every parsed point whose filename starts with `name`."""
        return {key: point for key, point in sorted(self.data.items()) if key.startswith(name)}

    def get_detector(self, name: str) -> np.ndarray:
        """Return one point's detector image, transposed to (row, column)."""
        camera_matrix = self.data[name]["Data"]["Detector"].T
        return camera_matrix

    def get_value(self, name: str, key: str) -> float:
        """Return one point's `key` value, e.g. '_2theta', searching every section of the file."""
        for section in self.data[name].values():
            if isinstance(section, dict) and key in section:
                return section[key]
        raise KeyError(f"{key!r} not found in {name!r}")

    def _attach_hkl_readout(self, ax: Axes, image: AxesImage, detector: np.ndarray, prefix: str = "") -> None:
        """Make the toolbar read-out of `ax` report hkl and counts for the pixel under the cursor."""
        rows, cols = detector.shape

        def format_coord(x: float, y: float) -> str:
            # imshow data coords are (column, row); round to the pixel under the cursor.
            j, i = int(round(x)), int(round(y))
            if not (0 <= i < rows and 0 <= j < cols):
                return f"row={y:.1f} col={x:.1f}"
            hkl = detector["hkl"][i, j]
            return f"{(j, i)}, ({hkl[0]:+.2f}, {hkl[1]:+.2f}, {hkl[2]:+.2f}), {detector['counts'][i, j]:+.2f}"

        ax.format_coord = format_coord
        # The toolbar appends the topmost mouseover artist's own bracketed value after format_coord's
        # string; dropping the image out of that set stops the count being reported a second time.
        image.set_mouseover(False)

    def show_detector(
        self, detector: np.ndarray, long_side_inches: float = 7.0, **imshow_kwargs: Any
    ) -> tuple[Figure, Axes]:
        """Plot the counts image with the cursor read-out reporting hkl and intensity per pixel."""
        rows, cols = detector.shape
        # Shape the canvas to the array, so square pixels fill it instead of leaving a sliver
        # in a default landscape figure. 1.8 in reserves room for the colorbar and labels.
        scale = long_side_inches / max(rows, cols)
        fig, ax = plt.subplots(figsize=(cols * scale + 1.8, rows * scale + 0.8), layout="constrained")
        image = ax.imshow(detector["counts"], origin="lower", cmap="viridis", **imshow_kwargs)
        fig.colorbar(image, ax=ax)
        self._attach_hkl_readout(ax, image, detector)
        plt.show()
        return fig, ax

    @staticmethod
    def _peak_title(label: str, peaks: list) -> str:
        """Caption a panel with its scan point and one line per peak found in it."""
        if not peaks:
            return f"{label}  no peak"
        lines = [label]
        for peak in peaks:
            row, col = peak["center"]
            h, k, el = peak["hkl"]
            lines.append(f"({row:.1f}, {col:.1f})  ({h:+.2f} {k:+.2f} {el:+.2f})  I={peak['counts']:.0f}")
        return "\n".join(lines)

    def show_scans(
        self,
        scan_name: str,
        ub_name: str,
        ncols: int = 6,
        panel_inches: float = 2.0,
        save: str | None = None,
        peak_kwargs: dict | None = None,
        **imshow_kwargs: Any,
    ) -> tuple[Figure, np.ndarray] | None:
        """Convert every point of `scan_name` to hkl and show them as a grid of plots in one window."""
        names = list(self.get_scans(scan_name))
        if not names:
            raise KeyError(f"no parsed points start with {scan_name!r}")
        detectors = [self.convert_camera_hkl(name, ub_name) for name in names]

        rows, cols = detectors[0].shape
        ncols = min(ncols, len(names))
        nrows = ceil(len(names) / ncols)
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(ncols * (panel_inches + 0.6) + 0.4, nrows * panel_inches * rows / cols + 0.6),
            layout="constrained",
            squeeze=False,
        )

        for ax, name, detector in zip(axes.ravel(), names, detectors):
            image = ax.imshow(detector["counts"], origin="lower", cmap="viridis", **imshow_kwargs)
            peaks = self.find_peaks(detector, **(peak_kwargs or {}))
            label = Path(name).stem.removeprefix(scan_name).lstrip("_")
            ax.set_title(self._peak_title(label, peaks), fontsize=5)
            for peak in peaks:
                row_slice, col_slice = peak["bounds"]
                # Outline the summed pixels themselves; the level sits between a masked pixel centre
                # and its unmasked neighbour, which is where the region's edge actually runs.
                ax.contour(
                    np.arange(col_slice.start, col_slice.stop),
                    np.arange(row_slice.start, row_slice.stop),
                    peak["mask"],
                    levels=[0.5],
                    colors="red",
                    linewidths=0.8,
                )
            ax.tick_params(labelsize=6)
            # Each point gets its own bar, so a weak point is not washed out by the brightest one.
            fig.colorbar(image, ax=ax).ax.tick_params(labelsize=6)
            self._attach_hkl_readout(ax, image, detector, prefix=f"{name}   ")
        for ax in axes.ravel()[len(names) :]:
            ax.set_axis_off()
        if save:
            # Write before show(): a blocking GUI backend only returns once the window is closed,
            # and the conversion above is too slow to repeat for a missing output directory.
            Path(save).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save, dpi=150)
        else:
            plt.show()
            return fig, axes

    def convert_camera_hkl(self, name: str, ub_name: str) -> np.ndarray:
        """Return one point as an (hkl, counts) record per pixel, using the UB stored as `ub_name`."""
        camera_matrix = self.get_detector(name)
        rows, cols = camera_matrix.shape
        dl = 116 / 512
        center = (255, 255)

        PIXEL_DTYPE = np.dtype([("hkl", np.float64, (3,)), ("counts", camera_matrix.dtype)])

        theta = -1 * self.get_value(name, "_2theta")
        phi_goni = -1 * self.get_value(name, "_phi")
        chi = -1 * self.get_value(name, "_chi")
        omega = -1 * self.get_value(name, "_omega")
        d = 409.9628
        ub = self.data["ub"][ub_name]

        # - for FIG.1 in ub pdf
        i, j = np.indices((rows, cols))
        xyz = np.stack([-(j - center[0]) * dl, (i - center[1]) * dl, np.zeros((rows, cols))], axis=-1)

        # transform_detector_coordinates, applied to the whole pixel grid at once.
        xyz = np.einsum("ij,...j->...i", rot_y(theta), xyz)
        # xyz = xyz @ rot_y(theta).T

        theta_rad = np.radians(theta)
        xyz[..., 2] += d * np.cos(theta_rad)
        xyz[..., 0] += d * np.sin(theta_rad)

        # q_l = k_i - k_f read straight off the pixel direction. This is the vector form of
        # get_hkl's (theta_pixel, phi_pixel) decomposition and agrees with it to 1e-15.
        k_neutron = 2 * np.pi / self.wavelength
        # norm converts xyz to unit vectors which is equivalent of converting to polar coordinates.
        n_hat = xyz / np.linalg.norm(xyz, axis=-1, keepdims=True)
        q_l = k_neutron * np.array([0.0, 0.0, 1.0]) - k_neutron * n_hat

        # r_mat and ub are constant across the detector, so invert them once rather than per pixel.
        r_mat = rot_y(omega) @ rot_z(chi) @ rot_y(phi_goni)
        transform = np.linalg.inv(ub) @ np.linalg.inv(r_mat) / (2 * np.pi)

        detector = np.empty((rows, cols), dtype=PIXEL_DTYPE)
        detector["hkl"] = q_l @ transform.T
        detector["counts"] = camera_matrix
        return detector

    def gaussian_filter_diff(
        self,
        detector: np.ndarray,
        sigma1: float | tuple[float, float],
        sigma2: float | tuple[float, float],
        **kwargs: Any,
    ) -> np.ndarray:
        """Return a copy of `detector` whose counts are the Gaussian blur at `sigma1` minus the one at `sigma2`."""
        # Counts come off the camera as integers, so the copy holds them as floats: the difference is
        # signed and well below one count away from a peak, which an integer field would round away.
        counts = detector["counts"].astype(np.float64)
        # Difference of Gaussians: with sigma1 < sigma2 the narrow blur keeps the peak while the wide one
        # estimates the background under it, so the difference leaves peaks standing on a flat zero.
        blur1 = ndimage.gaussian_filter(counts, sigma1, **kwargs)
        blur2 = ndimage.gaussian_filter(counts, sigma2, **kwargs)

        filtered = np.empty(detector.shape, dtype=[("hkl", np.float64, (3,)), ("counts", np.float64)])
        # hkl is where a pixel points, not what it measured; only the counts are filtered.
        filtered["hkl"] = detector["hkl"]
        filtered["counts"] = blur1 - blur2
        return filtered

    def find_peaks(
        self,
        detector: np.ndarray,
        sigma1: float = 8,
        sigma2: float = 16,
        threshold: float = 3.0,
        min_area: int = 9,
        min_counts: float = 1000,
    ) -> list[dict]:
        """Return every peak in `detector` as a list of dicts, strongest first, empty when it holds none."""
        counts = detector["counts"].astype(np.float64)
        dog = self.gaussian_filter_diff(detector, sigma1, sigma2)["counts"]
        # Cutting at a multiple of the image's own spread is self-regulating: real peaks inflate that
        # spread and so resolve into a few blobs, while a background-only frame fragments into dozens.
        labels, n = ndimage.label(dog > threshold * dog.std(), structure=np.ones((3, 3)))
        if n == 0:
            return []

        halo = int(round(2 * sigma2))  # the moat sits at r ~ sigma2, so look that far out
        found = []
        for lab, sl in enumerate(ndimage.find_objects(labels), start=1):
            if (labels[sl] == lab).sum() < min_area:
                continue
            pad = tuple(slice(max(0, s.start - halo), min(d, s.stop + halo)) for s, d in zip(sl, dog.shape))
            core = labels[pad] == lab
            # Restricting the ring to unlabelled pixels keeps a close neighbour's positive core from
            # filling in the moat and disqualifying a genuine peak.
            ring = ndimage.binary_dilation(core, np.ones((3, 3)), iterations=halo) & (labels[pad] == 0)
            surround = float(dog[pad][ring].mean()) if ring.any() else 0.0
            if surround > 0.0:  # a real peak is ringed by negative, noise is not
                continue
            weights = np.where(core, counts[pad], 0.0)
            total = float(weights.sum())
            if total < min_counts:
                continue
            # The raw counts set the centre, not the filtered image: the blur that found the blob also
            # drags its shape about, while the counts are what the reflection actually deposited.
            row, col = ndimage.center_of_mass(weights)
            row, col = row + pad[0].start, col + pad[1].start
            # The centre lands between pixels, so hkl is interpolated rather than read off the nearer one.
            hkl = np.array(
                [ndimage.map_coordinates(detector["hkl"][..., k], [[row], [col]], order=1)[0] for k in range(3)]
            )
            found.append(
                {
                    "center": (row, col),
                    "hkl": hkl,
                    "counts": total,
                    "bounds": sl,
                    # The summed region is ragged, so it is carried as a mask rather than inferred from
                    # `bounds`: the box around it holds about half again as many pixels.
                    "mask": labels[sl] == lab,
                    "surround": surround,
                }
            )
        return sorted(found, key=lambda peak: -peak["counts"])

    # def convert_camera_hkl(self, name, ub_name):
    #     camera_matrix = self.get_detector(name)
    #     rows, cols = camera_matrix.shape
    #     dl = 116 / 512
    #     center = (255, 255)

    #     PIXEL_DTYPE = np.dtype([("hkl", np.float64, (3,)), ("counts", camera_matrix.dtype)])

    #     theta = -1 * self.get_value(name, "_2theta")
    #     phi_goni = -1 * self.get_value(name, "_phi")
    #     chi = -1 * self.get_value(name, "_chi")
    #     omega = -1 * self.get_value(name, "_omega")
    #     d = 409.9628
    #     ub = self.data["ub"][ub_name]
    #     detector = np.empty((rows, cols), dtype=PIXEL_DTYPE)
    #     for i in range(rows):
    #         for j in range(cols):
    #             # - for FIG.1 in ub pdf
    #             x = -(j - center[0]) * dl
    #             y = (i - center[1]) * dl
    #             z = 0
    #             xyz = np.array([x, y, z])
    #             detector_xyz = self.transform_detector_coordinates(xyz,theta, d)

    #             theta_pixel, phi_pixel = self.get_polar_coordinate(detector_xyz)

    #             hkl = self.get_hkl(self.wavelength, theta_pixel, phi_pixel, phi_goni, chi, omega, ub)
    #             detector[i, j] = (hkl, camera_matrix[i, j])
    #     return detector
    # def transform_detector_coordinates(self,  xyz, theta, d):
    #     """
    #     DEMAND Convention -> Mantid: theta(or 2theta), phi, omega, chi all need * -1.
    #     """

    #     # Rotate detector by theta (2theta)
    #     xyz = np.asarray(xyz)
    #     new_xyz = rot_y(theta) @ xyz.T

    #     # translate along theta(2theta by d), y is height doesn't change.
    #     x, y, z = new_xyz
    #     theta_rad = np.radians(theta)
    #     z += d * np.cos(theta_rad)
    #     x += d * np.sin(theta_rad)
    #     return np.array([x, y, z])

    # def get_polar_coordinate(self, xyz):
    #     """Get theta and phi"""
    #     x, y, z = xyz
    #     cos_theta = np.dot(xyz, [0,0,1]) / np.linalg.norm(xyz)
    #     theta = np.degrees(np.arccos(cos_theta))
    #     theta *= -1 if x < 0 else 1
    #     phi = np.degrees(np.arctan(y/x))
    #     return theta, phi

    # def get_hkl(self, wavelength, theta, phi_pixel, phi_goni,chi, omega, ub):
    #     """
    #     theta, phi is calculated, so already satisfy mantid convention. chi, omega is read from data, need to *= -1.
    #     """
    #     r_mat = rot_y(omega) @ rot_z(chi) @ rot_y(phi_goni)
    #     k_neutron = 2 * np.pi / wavelength
    #     theta_rad = np.radians(theta)
    #     phi_pixel_rad = np.radians(phi_pixel)
    #     q_l = [-k_neutron * np.sin(theta_rad) * np.cos(phi_pixel_rad),
    #         -k_neutron * np.sin(theta_rad) * np.sin(phi_pixel_rad),
    #         k_neutron - k_neutron * np.cos(theta_rad)]
    #     hkl = np.linalg.inv(ub) @ np.linalg.inv(r_mat) @ q_l / (2 * np.pi)
    #     return hkl
