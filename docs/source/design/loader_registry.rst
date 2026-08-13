LoaderRegistry and Loaders
==========================

Overview
--------

The **LoaderRegistry** is a singleton that manages the registration and lifecycle of data loaders in TAVI. Loaders are responsible for parsing and loading scan data from various file formats.

.. _loader-architecture:

Architecture
------------

**LoaderRegistry** (``tavi.library.storage.loader.loader_registry``)
    A singleton registry (``@Singleton`` from
    ``neutrons_standard.decorators.singleton``) that:

    - Registers the built-in loaders (``ORNLSpiceLoader``, ``DefaultLoader``) on
      construction, keyed by each loader's ``get_scan_type()``
    - Propagates the filestore to all registered loaders
    - Provides access to loaders via ``get_loaders()`` and ``get_loader(key)``

    The constructor takes **no arguments** and starts with ``filestore = None``.
    The filestore is injected afterwards by ``tavi.library.init()``:

    .. code-block:: python

        filestore = LocalFileStore()
        loader_registry = LoaderRegistry()
        loader_registry.set_filestore(filestore)
        RawScanLoadController(filestore)

**AbstractLoader** (``tavi.library.storage.loader.interface.base``)
    Base class for all loaders. Implements:

    - ``__init__(filestore)``: Store the filestore reference
    - ``set_filestore(filestore)``: Update the filestore reference
    - ``generate_uuid(file_path) -> UUID``: MD5 of the file's text content
    - Inherits abstract methods from ``LoaderInterface``

**LoaderInterface** (``tavi.library.storage.loader.interface.loader_interface``)
    Abstract interface defining loader contract:

    - ``load(path: str) -> Scan``: Load complete scan data
    - ``get_scan_type() -> RawScanType``: Registry key / classification result
    - ``get_score(path: str) -> float``: Assess loader suitability for a file
    - ``parse_metadata(path: str) -> ScanMetadata``: Extract metadata
    - ``parse_tavi_metadata(path: str) -> TaviMetadata``: Extract TAVI-specific
      metadata (default axes, friendly name/path, normalization)
    - ``parse_scan_values(path: str) -> ScanData``: Extract scan values
    - ``parse_external_metadata(path: str) -> dict``: Extract additional metadata
      from sibling files (for SPICE, the ``UBConf`` entry)
    - ``adapt_scan_data(meta, tavi_meta, values) -> RawScan``: Combine parsed data
      into a scan object

Usage
-----

Registering a Loader
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.data.enum.raw_scan_type import RawScanType
    from tavi.library.storage.loader.loader_registry import LoaderRegistry
    from my_loaders import MyCustomLoader

    registry = LoaderRegistry()
    registry.register(RawScanType.MyFormat, MyCustomLoader(filestore))

The key is a :class:`tavi.library.data.enum.raw_scan_type.RawScanType` member —
the same value ``RawScanClassifier`` returns and ``RawScanLoadController`` looks
up. ``_register_loader(loader)`` is the internal shorthand that derives the key
from ``loader.get_scan_type()``.

Retrieving Loaders
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    loaders = registry.get_loaders()          # all registered loaders
    loader = registry.get_loader(classification)  # one loader by RawScanType

Both refresh the filestore on every registered loader before returning.
``get_loader`` raises ``RuntimeError(f"No loader for key: {key}")`` if the key
was never registered.

Updating Filestore
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    registry.set_filestore(new_filestore)  # updates all loaders automatically

Implementing a Custom Loader
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.data.enum.raw_scan_type import RawScanType
    from tavi.library.storage.loader.interface.base import AbstractLoader

    class MyLoader(AbstractLoader):
        def get_scan_type(self) -> RawScanType:
            return RawScanType.MyFormat

        def load(self, path: str) -> Scan:
            # Implementation
            pass

        def get_score(self, path: str) -> float:
            # Return score indicating suitability for file
            pass

        # Implement the remaining abstract methods...

A full walkthrough, including the enum and registry edits, is in
:doc:`../guides/loader_setup`.

Built-in Loaders
----------------

``ORNLSpiceLoader`` (``RawScanType.ORNLSpice``)
    Loads ORNL SPICE ``.dat`` files. Scores files with a
    :doc:`../guides/rule_based_classifier` over ``ORNLSpiceRuleSet``.

``DefaultLoader`` (``RawScanType.NONE``)
    Fallback. Always scores ``0`` and raises
    ``RuntimeError(f"No suitable loader found for file at: {path}")`` if anything
    actually tries to load through it. Its role is to give the classifier a
    well-defined "nothing matched" answer rather than to load anything.

Key Features
------------

- **Singleton Pattern**: Only one registry instance exists application-wide
- **Filestore Synchronization**: ``get_loaders``/``get_loader``/``set_filestore``
  push the current filestore to every loader
- **Extensible Design**: Easy to add new loaders by implementing ``AbstractLoader``
- **Score-based Selection**: Loaders indicate suitability via ``get_score()``;
  ``RawScanClassifier`` picks the highest scorer

See Also
--------

- :doc:`taviclasses` - Core TAVI class structure
- :doc:`raw_scan_load_controller` - The controller that drives classification and loading
- :doc:`filestore` - The file storage interface loaders read through
