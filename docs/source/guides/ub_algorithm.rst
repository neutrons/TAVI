.. _ub_algorithm:

++++++++++++
UB Algorithm
++++++++++++

``UBAlgorithm`` implements the orientation-matrix (UB) determination used by
:ref:`TAS <triple_axis>`. It builds the U matrix from indexed peaks and the
sample's B matrix, supporting one-, two-, three-, and many-peak determinations.


Overview
========

Source module:

.. code-block:: text

   tavi.library.ubalgorithm.ub

Primary class:

.. autoclass:: tavi.library.ubalgorithm.ub.UBAlgorithm
   :members:
   :undoc-members:


Method Notes
============

``find_u_from_one_peak_and_scattering_plane``
---------------------------------------------

* Determines U from a single peak plus a defined scattering plane. The second
  in-plane reference is chosen to be as orthogonal as possible to the peak.

``find_u_from_two_peaks``
-------------------------

* Determines U from two non-collinear peaks, assuming the cross product defines
  the out-of-plane direction.

``find_ub_from_three_peaks`` / ``find_ub_from_multiple_peaks``
--------------------------------------------------------------

* Determine the full UB from three peaks, or a least-squares UB from many peaks.

``plane_normal_from_two_peaks``
-------------------------------

* Returns the scattering-plane normal implied by two peaks.


Coordinate Conventions
======================

The library computes UB in the Mantid frame. Use
:func:`tavi.library.experiment.utilities.mantid_to_spice` (and its inverse
``spice_to_mantid``) to convert to/from the SPICE frame. Both accept a
``version`` argument (``"new"`` by default, ``"old"`` for the legacy transform).


Minimal Example
===============

.. code-block:: python

   from tavi.library.experiment.utilities import mantid_to_spice

   # tas.ub returns the UB matrix in the Mantid frame.
   ub_mantid = tas.ub(peaks=(peaks[0], peaks[1]))
   ub_spice = mantid_to_spice(ub_mantid)
