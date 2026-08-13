.. _scan_browser:

++++++++++++
Scan Browser
++++++++++++

``browse_scans`` draws a grid of scans -- one subplot per scan number -- with
optional fits, per-component curves, and resolution bars. It is the rendering
half of :ref:`TAS.browse <triple_axis>`; ``TAS`` computes the resolution and
forwards everything else straight through.


Overview
========

Source module:

.. code-block:: text

   tavi.library.plot.browser

.. autofunction:: tavi.library.plot.browser.browse_scans

The resolution-ellipse plots live separately in
``tavi.library.plot.plot_ellipse.PlotResolution``, which now holds only the
ellipse rendering:

.. autoclass:: tavi.library.plot.plot_ellipse.PlotResolution
   :members:
   :undoc-members:


Layout
======

The grid is at most three columns wide (``ncols = min(3, len(scan_list))``) with
as many rows as needed; leftover axes in the last row are hidden. Each subplot
is titled ``"<scan number>, <hkl>"``.

The x-axis defaults to ``scan.metadata.def_x`` and the y-axis to
``scan.metadata.def_y`` unless ``def_x`` / ``def_y`` are given. Column names
starting with a digit are prefixed with an underscore to match the loader's
attribute-safe keys (for example ``2theta`` becomes ``_2theta``). Counts are
drawn with :math:`\sqrt{y}` error bars.

``xlim`` and ``ylim`` accept either form:

* a **scalar** -- scales the axis symmetrically about the data midpoint, so
  ``xlim=1.1`` widens the range to 1.1x the data span;
* a **[min, max] list** -- sets the range directly.


Fits and Components
===================

With ``show_fits=True`` each scan is fit with ``model_dict`` through the
selected :ref:`FitPackage <fit_package>`, and the composite curve is overlaid.
The legend label carries the fitted FWHM of every peak, ``fit=<fwhm>, <fwhm>``.

``show_components=True`` additionally draws each component of a multi-component
model as a dashed line, labelled by its prefix. A component that has a ``fwhm``
parameter is labelled ``<prefix>=<fwhm>``; because the widths then appear on the
individual components, the composite curve's label collapses to plain ``fit``.
Single-component fits are unaffected -- the component overlay only kicks in when
``len(fit_result.components) > 1``.


.. _resolution_bar:

Resolution Bars
===============

``resolution_bars`` is the value returned by
:ref:`TAS.resolution_bar <triple_axis>`, i.e. a 2-tuple
``(coherent_fwhm_per_scan, res_4d_per_scan)``. Passing it switches on three
behaviours:

**A horizontal bar per fitted peak.** Each bar is centred on the peak's fitted
centre, has a total width equal to the coherent FWHM, and is drawn at the peak's
half-maximum. The half-maximum is measured from the background: components
without a fitted ``center`` (a linear background, say) are evaluated at the peak
centre and added underneath, so ``bar_y = background + height / 2``. Only the
first bar of a subplot is labelled, ``res=<fwhm>``, with one entry per peak.

**Multiple peaks per scan.** The coherent FWHM for a scan is a tuple aligned
with ``fit_result.peaks``, so a scan fit with two Gaussians gets two bars of
independent widths. Bar ``idx`` uses width ``coh_list[idx]``.

**An x-axis in** :math:`\Delta q`. When the default x motor is an angle
(``s1``, ``s2`` or ``omega``) the abscissa is replaced by
``experiment.get_delta_q(...)`` and relabelled ``del_q(<motor>)``, because the
bar width is a :math:`Q` width and would be meaningless against degrees. Scans
whose default x is already a non-angle motor (an energy scan, for instance) keep
their own axis and the bar is drawn in those units.

``show_fits`` must be ``True`` -- the bar is positioned from the fit, so
``browse_scans`` raises ``ValueError`` if resolution bars are requested without
fits.


Return Value
============

``browse_scans`` always returns a 3-tuple. Without resolution bars every element
is ``None``:

.. code-block:: python

   hkls, fit_results, res_mat_4d = browse_scans(...)                       # (None, None, None)
   hkls, fit_results, res_mat_4d = browse_scans(..., resolution_bars=bars) # populated

* ``hkls`` -- ``(h, k, l)`` per scan, in ``scan_list`` order,
* ``fit_results`` -- the ``FitResult`` per scan (see :ref:`fit_package`),
* ``res_mat_4d`` -- the 4D resolution matrices and ``r0`` values passed in, one
  tuple of ``(res_mat, r0)`` per fitted peak per scan.

Together these are everything needed for a downstream intensity export:
positions, integrated intensities, and the resolution volume to normalise by.
