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

``calculate_resolution``
------------------------

* Builds the ``Resolution`` object (default Cooper-Nathans) used by downstream
  resolution-ellipse calculations.

``browse`` / ``browse_resolution_ellipse`` / ``resolution_bar``
---------------------------------------------------------------

* Plot a grid of scans (optionally with fits and resolution bars) and visualize
  per-scan resolution ellipses. ``fit_package`` and ``model_dict`` control how
  each scan is fit.


Minimal Example
===============

.. code-block:: python

   from tavi.library.Instrument.instrument import Instrument, InstrumentCatalog
   from tavi.library.experiment.experiment import Experiment
   from tavi.library.experiment.enum import FixedEnergyMode
   from tavi.library.tas.triple_axis import TAS

   experiment = Experiment(mode=FixedEnergyMode.FIX_Ef, ei_or_ef=14.503)
   sample = experiment.create_sample("/path/to/UBConf/UBmatrix.dat")
   tas = TAS(Instrument(instrument_catalog=InstrumentCatalog.HB1A4C), sample, experiment)

   # Refine the UB matrix from two indexed peaks.
   ub = tas.ub(peaks=(peaks[0], peaks[1]))
