"""Data browser."""

from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.fit import FitPackage
from tavi.library.fit.fit import ModelName


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
    from tavi.library.fit import Fit

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
    print(bars)
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
        if show_resolution_bar and def_x in ["s1", "s2", "omega"]:
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

            # When components are shown separately, each component carries its own
            # fwhm in its label, so the composite fit label is kept as plain "fit".
            show_component_labels = show_components and len(fit_result.components) > 1
            if show_component_labels:
                fit_label = "fit"
            else:
                fwhm_label = ", ".join(f"{peak.values['fwhm']:.3g}" for peak in peaks)
                fit_label = f"fit={fwhm_label}"
            ax.plot(x_fine, fit_result.raw.eval(x=x_fine), label=fit_label)

            # Optionally overlay each component (peak, linear background, ...) separately.
            if show_component_labels:
                for prefix, y_comp in fit_result.raw.eval_components(x=x_fine).items():
                    label = prefix.rstrip("_") or "component"
                    component = fit_result.components.get(prefix) or fit_result.components.get("")
                    if component is not None and "fwhm" in component.values:
                        label = f"{label}={component.values['fwhm']:.3g}"
                    ax.plot(x_fine, y_comp, "--", lw=1, label=label)

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

                coh_list = np.atleast_1d(coh)
                if show_component_labels:
                    sep = "\n" if len(coh) > 1 else ", "
                    reso_label = "res=" + sep.join(f"{c:.3g}" for c in coh_list)
                else:
                    reso_label = "res=" + ",".join(f"{c:.3g}" for c in coh_list)
                for idx, peak in enumerate(peaks):
                    center = peak.values["center"]
                    components = fit_result.raw.eval_components(x=center)
                    background = sum(float(value) for key, value in components.items() if is_background(key))
                    bar_y = background + peak.values["height"] / 2
                    ax.errorbar(
                        center,
                        bar_y,
                        xerr=coh_list[idx] / 2,
                        color="red",
                        capsize=4,
                        lw=2,
                        label=reso_label if idx == 0 else None,
                    )
            ax.legend(fontsize=8, loc="upper right")

        # if we can't parse title for nominal hkl position then we use a Gaussian to try to find the center of the peak.
        try:
            hkl = experiment.get_hkl(dict(scan_num=num), use_title=True)
        except ValueError:
            hkl = experiment.get_hkl(
                dict(scan_num=num), use_title=False, model_dict=[(ModelName.Gaussian, dict(guess=True))]
            )
        hkls.append(hkl)
        ax.set_title(f"{num}, {hkl}")
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
    return None, None, None
