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

``ORNLSpiceLoader.load(file_path)`` performs these steps:

1. Generate the scan's uuid via ``generate_uuid`` (MD5 of the file's text).
2. Parse numeric scan arrays via ``parse_scan_values``.
3. Parse header metadata via ``parse_metadata``.
4. Parse TAVI metadata via ``parse_tavi_metadata``.
5. Build provenance via ``create_provenance``.
6. Read the ``ubconf`` field out of the parsed metadata and resolve it with
   ``parse_external_metadata(file_path, ub_name)``.
7. Merge the external UB data into ``meta.data`` and return a ``RawScan`` via
   ``adapt_scan_data``.

Note that this loader's ``parse_external_metadata`` takes a second argument
(``ub_name``) and its ``adapt_scan_data`` takes ``uuid``, ``values``, ``meta``,
``tavi_meta`` and ``prov`` — both are wider than the corresponding
``LoaderInterface`` signatures. ``load()`` is the only caller, so the registry
never sees the difference.


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
* Groups all fields under a single ``"ORNL Metadata"`` category (see
  :doc:`../design/frontend/data_file_view`) so ``.``-delimited access
  (``scan.metadata.ubconf``) keeps working, while the metadata widget renders
  them under one ``"ORNL Metadata"`` tab.

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
