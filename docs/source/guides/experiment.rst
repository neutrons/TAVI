.. _experiment:

++++++++++
Experiment
++++++++++

``Experiment`` is the data-access layer that sits between a loader (for example
:ref:`ORNLSpiceLoader <ornl_spice_loader>`) and the analysis classes. It loads
scans from disk, looks them up by scan number, and derives per-scan quantities
such as peak centers, ``(h, k, l)``, and momentum transfer.


Overview
========

Source module:

.. code-block:: text

   tavi.library.experiment.experiment

Primary class:

.. autoclass:: tavi.library.experiment.experiment.Experiment
   :members:
   :undoc-members:


Fixed-Energy Mode
=================

An ``Experiment`` is created in a fixed-energy mode (fixed :math:`E_i` or fixed
:math:`E_f`) with the corresponding energy value. That mode determines how
``get_ei_ef`` maps a fitted energy transfer back to incident/scattered energies,
and which energy is held constant when peaks are located.


Method Notes
============

``load_file`` / ``load_folder``
-------------------------------

* Ingest a single scan file or every scan in a directory into the experiment.

``get_peak_center``
-------------------

* Fits the scan(s) and returns a ``list[DataPoint]`` -- one per fitted peak
  center. Inelastic scans may yield more than one ``DataPoint``.
* Takes a ``FitPackage`` and a ``model_dict`` (see :ref:`fit_package`).

``get_closest_to_center_data_point``
------------------------------------

* Like ``get_peak_center``, but returns the *measured* data point nearest each
  fitted centre rather than the fitted centre itself, so the motor angles are
  real ones the goniometer can turn into a rotation matrix. Used by
  :ref:`TAS.resolution_bar <triple_axis>`.
* The experiment's ``mode`` and its corresponding fixed energy are forwarded to
  the loader, so each returned ``DataPoint`` carries a consistent
  :math:`(E_i, E_f)` pair.

``get_hkl`` / ``get_delta_q`` / ``get_two_theta`` / ``get_psi``
---------------------------------------------------------------

* Derive reciprocal-space and angular quantities for a located scan.
* ``get_delta_q`` also supplies the abscissa used by the
  :ref:`scan browser <scan_browser>` whenever a resolution bar is drawn against
  an angular motor.

``create_sample``
-----------------

* Builds a ``Sample`` (oriented lattice + UB) from a UB configuration file.


Minimal Example
===============

.. code-block:: python

   from tavi.library.experiment.experiment import Experiment
   from tavi.library.experiment.enum import FixedEnergyMode
   from tavi.library.fit.fit import FitPackage, ModelName

   experiment = Experiment(mode=FixedEnergyMode.FIX_Ef, fixed_energy=14.503)
   experiment.load_folder("/path/to/exp1091/Datafiles")

   # One DataPoint per fitted peak; flatten across scans into a single list.
   scan_list = [531, 532, 533]
   peaks = [
       dp
       for i in scan_list
       for dp in experiment.get_peak_center(
           {"scan_num": i}, FitPackage.lmfit, [(ModelName.Gaussian, dict(guess=True))]
       )
   ]
   print(peaks[0].hkl, peaks[0].ei, peaks[0].ef)
