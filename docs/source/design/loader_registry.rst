LoaderRegistry and Loaders
==========================

Overview
--------

The **LoaderRegistry** is a singleton that manages the registration and lifecycle of data loaders in TAVI. Loaders are responsible for parsing and loading scan data from various file formats.

.. _loader-architecture:

Architecture
------------

**LoaderRegistry** (``tavi.library.storage.loader.loader_registry``)
    A singleton registry that:

    - Maintains a registry of all available loaders
    - Propagates the filestore to all registered loaders
    - Provides access to loaders via ``get_loaders()``

**AbstractLoader** (``tavi.library.storage.loader.interface.base``)
    Base class for all loaders. Implements:

    - ``set_filestore(filestore)``: Update the filestore reference
    - Inherits abstract methods from ``LoaderInterface``

**LoaderInterface** (``tavi.library.storage.loader.interface.loader_interface``)
    Abstract interface defining loader contract:

    - ``load(path: str) -> Scan``: Load complete scan data
    - ``get_score(path: str) -> int``: Assess loader suitability for a file
    - ``parse_metadata(path: str) -> ScanMetadata``: Extract metadata
    - ``parse_ub(path: str) -> ScanUb``: Extract UB matrix data
    - ``parse_scan_values(path: str) -> ScanValues``: Extract scan values
    - ``parse_external_metadata(path: str) -> dict``: Extract additional metadata
    - ``adapt_scan_data(...) -> Scan``: Combine parsed data into Scan object

Usage
-----

Registering a Loader
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.storage.loader.loader_registry import LoaderRegistry
    from my_loaders import MyCustomLoader

    registry = LoaderRegistry(filestore)
    custom_loader = MyCustomLoader(filestore)
    registry.register("MyLoader", custom_loader)

Retrieving Loaders
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    loaders = registry.get_loaders()  # Returns all registered loaders
    # Filestore is automatically refreshed for all loaders

Updating Filestore
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    registry.set_filestore(new_filestore)  # Updates all loaders automatically

Implementing a Custom Loader
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.storage.loader.interface.base import AbstractLoader

    class MyLoader(AbstractLoader):
        def load(self, path: str) -> Scan:
            # Implementation
            pass

        def get_score(self, path: str) -> int:
            # Return score indicating suitability for file
            pass

        # Implement other abstract methods...

Key Features
------------

- **Singleton Pattern**: Only one registry instance exists application-wide
- **Filestore Synchronization**: Changes to filestore propagate to all loaders
- **Extensible Design**: Easy to add new loaders by implementing ``AbstractLoader``
- **Score-based Selection**: Loaders can indicate suitability via ``get_score()``

See Also
--------

- :doc:`taviclasses` - Core TAVI class structure
- File storage interface documentation
