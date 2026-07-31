.. _ornl_spice_loader:

+++++++++++++++
ORNLSpiceLoader
+++++++++++++++

``ORNLSpiceLoader`` loads ORNL SPICE scan files (``.dat``) into TAVI scan models.
It parses:

* scan values from the tabular section,
* scan metadata from file headers,
* TAVI-specific metadata used by the application layer,
* optional UB configuration data from sibling ``UBConf`` files.


Overview
========

Source module:

.. code-block:: text

   tavi.library.storage.loader.ornl_spice_loader

Primary class:

.. autoclass:: tavi.library.storage.loader.ornl_spice_loader.ORNLSpiceLoader
   :members:
   :undoc-members:


Expected Data Layout
====================

The loader expects a structure similar to:

.. code-block:: text

   exp815/
     Datafiles/
       HB1_exp0815_scan0003.dat
     UBConf/
       UB13Jun2019_41635PM.ini

The ``ubconf`` entry in the ``.dat`` header is used to resolve the UB file under
the sibling ``UBConf`` directory.


Typical Loading Flow
====================

``ORNLSpiceLoader.load(path)`` performs these steps:

1. Parse numeric scan arrays via ``parse_scan_values``.
2. Parse header metadata via ``parse_metadata``.
3. Parse TAVI metadata via ``parse_tavi_metadata``.
4. Build provenance via ``create_provenance``.
5. Resolve and parse external UB metadata via ``parse_external_metadata``.
6. Merge external UB data into metadata and return a ``RawScan``.


Method Notes
============

``parse_scan_values``
---------------------

* Converts SPICE column names into valid attribute-like keys.
* Handles empty/invalid measurement sections by returning empty arrays.

``parse_metadata``
------------------

* Collects key/value header fields before ``# col_headers =``.
* Captures terminal status entries (for example scan completed/stopped).

``parse_external_metadata``
---------------------------

* Supports INI-like UB files and legacy XML UB content.
* Returns a plain dictionary that can be merged into scan metadata.

``get_data_point_closest_to_center``
------------------------------------

* Fits the scan with ``model_dict`` and returns a ``list[DataPoint]`` -- one per
  fitted centre -- for the measured row nearest each centre. Multi-peak models
  therefore yield several data points from a single scan.
* Takes ``mode`` (``FixedEnergyMode``) and ``fixed_energy`` so :math:`E_i` and
  :math:`E_f` can be reconstructed per point. ``Experiment`` supplies both from
  its own mode.
* When the default x column is an angle (``s1``, ``s2``, ``omega``) the fit runs
  against :math:`\Delta q` from ``get_delta_q`` instead of the raw motor, so the
  centre is located in the same units the resolution bar is drawn in. Other
  default-x columns are fit against their own values.
* The ``MotorAngles`` attached to each ``DataPoint`` uses the goniometer's key
  names -- ``two_theta`` and ``omega`` rather than the SPICE column names ``s2``
  and ``s1``, plus ``sgl`` and ``sgu`` -- so ``Goniometer.r_mat`` can consume
  them directly when projecting resolution into an ``hkle`` frame. See
  :ref:`resolution_rotation_matrix`.


Minimal Example
===============

.. code-block:: python

   from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader
   from tavi.library.storage.local_file_store import LocalFileStore

   loader = ORNLSpiceLoader(LocalFileStore())
   scan = loader.load("test_data/exp815/Datafiles/HB1_exp0815_scan0003.dat")

   print(scan.uuid)
   print(scan.metadata.ubconf)
   print(scan.data.h[:3])
