.. _triple_axis:

+++++++++++++++++++
TAS (Triple-Axis)
+++++++++++++++++++

``TAS`` is the top-level entry point of the TAVI library. It ties together an
``Instrument``, a ``Sample``, and an :ref:`Experiment <experiment>`, and exposes
the high-level operations: UB-matrix determination, resolution calculation, and
scan browsing.


Overview
========

Source module:

.. code-block:: text

   tavi.library.tas.triple_axis

Primary class:

.. autoclass:: tavi.library.tas.triple_axis.TAS
   :members:
   :undoc-members:

UB convention:

.. autoclass:: tavi.library.tas.triple_axis.UBConvention
   :members:
   :undoc-members:


Method Notes
============

``ub``
------

* Determines the UB matrix from 1, 2, 3, or more peaks. With a single peak a
  ``scattering_plane`` must be supplied. Delegates to
  :ref:`UBAlgorithm <ub_algorithm>`.

``resolution_bar``
------------------

* Computes the coherent FWHM for every fitted peak of every scan in
  ``scan_list``, in the frame given by ``resolution_frame`` and along
  ``projection_axis``. See :ref:`resolution` for what those two arguments mean.
* Builds its own ``Resolution`` object from ``resolution_model`` (default
  Cooper-Nathans), so no separate setup call is needed.
* Peak centres come from
  ``Experiment.get_closest_to_center_data_point``, which fits each scan with
  ``model_dict`` and returns the measured data point nearest each fitted centre.
  Only ``ORNLSpiceLoader`` data is supported; anything else raises
  ``ValueError("Data format not implemented yet.")``.
* Returns ``(resolution_bar, res_4ds)``. Both are lists with one entry per scan,
  and each entry is a tuple with one item per fitted peak in that scan:

  .. code-block:: python

     resolution_bar = [(fwhm_peak0, fwhm_peak1), (fwhm_peak0,), ...]
     res_4ds        = [((res_mat, r0), (res_mat, r0)), ((res_mat, r0),), ...]

  The per-scan grouping is what lets the browser draw one bar per peak; a flat
  list would lose the scan boundaries.
* For ``hkle`` frames the incoherent FWHM is computed as well and printed as
  ``incoherent FWHM = [...]``. It is not returned or plotted yet.

``browse``
----------

* Plots the scans through :ref:`browse_scans <scan_browser>`, passing
  ``show_fits``, ``model_dict``, ``show_components``, ``def_x``/``def_y`` and
  ``xlim``/``ylim`` straight through.
* ``show_resolution_bar=True`` first calls ``resolution_bar`` with the same
  ``resolution_frame``, ``projection_axis`` and ``model_dict``, then hands the
  result to the browser. In that case ``browse`` returns
  ``(hkls, fit_results, res_mat_4d)`` for intensity export; otherwise it returns
  ``(None, None, None)``. It always returns a 3-tuple, so unpacking is safe
  either way — but truth-testing the result is not, since a 3-tuple of ``None``
  is itself truthy.
* ``show_fits`` must stay ``True`` when the resolution bar is on -- the bar is
  positioned from the fitted centre.

``browse_resolution_ellipse``
-----------------------------

* Plots the :math:`(Q_\parallel, Q_\perp)` resolution ellipse per scan, always in
  the local :math:`Q` frame, annotated with the coherent FWHM along axis 0
  (parallel) and axis 1 (perpendicular). Rendering is handled by
  ``PlotResolution.plot_resolution_ellipse``.


.. _tas_resolution_frame_arguments:

Resolution Arguments
====================

``resolution_frame`` and ``projection_axis`` appear on both ``browse`` and
``resolution_bar`` and are forwarded unchanged to the resolution stack.

``resolution_frame``
   ``"local"`` (the default) for the frame attached to the scattering triangle,
   or a 4-tuple such as ``((1, 0, 0), (0, 1, 0), (0, 0, 1), "e")`` for an
   ``hkle`` frame. A value that is neither raises
   ``ValueError("Resolution frame is not defined properly.")``.

``projection_axis``
   Index into the chosen frame: ``0``/``1``/``2`` for the three momentum axes
   and ``3`` for energy. In the local frame ``0`` is along :math:`\vec Q` and
   ``1`` is the in-plane perpendicular; in an ``hkle`` frame the indices follow
   the order the projection vectors were listed in.

Full details, including the axis conventions and how the sample rotation matrix
is obtained, are in :ref:`resolution`.

.. note::

   ``TAS.calculate_resolution`` has been removed. ``resolution_bar`` and
   ``browse_resolution_ellipse`` now construct their own ``Resolution`` with the
   frame they need, so ``tas.resolution`` is populated as a side effect of
   calling them rather than by a separate setup step. Code that called
   ``tas.calculate_resolution(...)`` first can simply drop the call.


Minimal Example
===============

.. code-block:: python

   from tavi.library.Instrument.instrument import Instrument, InstrumentCatalog
   from tavi.library.experiment.experiment import Experiment
   from tavi.library.experiment.enum import FixedEnergyMode
   from tavi.library.tas.triple_axis import TAS

   experiment = Experiment(mode=FixedEnergyMode.FIX_Ef, fixed_energy=14.503)
   sample = experiment.create_sample("/path/to/UBConf/UBmatrix.dat")
   tas = TAS(Instrument(instrument_catalog=InstrumentCatalog.HB1A4C), sample, experiment)

   # Refine the UB matrix from two indexed peaks.
   ub = tas.ub(peaks=(peaks[0], peaks[1]))


Browsing With a Resolution Bar
==============================

Peak plus linear background, energy resolution in a reordered ``hkle`` frame:

.. code-block:: python

   from tavi.library.fit import ModelName

   hkls, fit_results, res_4d = tas.browse(
       scan_list,
       show_fits=True,
       model_dict=[
           (ModelName.Gaussian, dict(prefix="g", guess=True)),
           (ModelName.Linear, dict(prefix="l", guess=True)),
       ],
       show_resolution_bar=True,
       show_components=False,
       resolution_frame=((1, 0, 0), (0, 0, 1), (0, 1, 0), "e"),
       projection_axis=3,
   )

The same call in the local :math:`Q` frame, taking the width along
:math:`\vec Q`, needs neither argument -- ``resolution_frame`` defaults to
``"local"`` and ``projection_axis`` to ``0``:

.. code-block:: python

   hkls, fit_results, res_4d = tas.browse(
       scan_list,
       model_dict=[(ModelName.Gaussian, dict(guess=True))],
       show_resolution_bar=True,
   )
