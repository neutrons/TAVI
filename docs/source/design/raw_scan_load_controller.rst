RawScanLoadController
=====================

.. class:: RawScanLoadController(filestore: FileStoreInterface)

   Orchestrates the loading of :class:`RawScan` objects from a file store.
   This controller centralizes loader selection, classification, and batch-loading
   workflows. Implemented as a singleton to ensure consistent loader registry
   and classifier usage across the application.

   :param filestore: File storage backend used to resolve file paths.
   :type filestore: FileStoreInterface

Overview
--------

The controller coordinates three main concerns:

- **Classification**: Determines the appropriate loader for a file via
  ``RawScanClassifier``.
- **Loader resolution**: Uses ``LoaderRegistry`` to map classifications to
  concrete loader implementations.
- **Batch orchestration**: Provides optimized loading for multiple files or folders.

Both ``LoaderRegistry`` and ``RawScanClassifier`` are constructed with no
arguments — they are singletons resolved at ``__init__`` time. ``Config`` and
``Singleton`` come from ``neutrons_standard``:

.. code-block:: python

    from neutrons_standard.config import Config
    from neutrons_standard.decorators.singleton import Singleton

    @Singleton
    class RawScanLoadController:
        def __init__(self, filestore: FileStoreInterface) -> None:
            self.filestore = filestore
            self.loader_registry = LoaderRegistry()
            self.classifier = RawScanClassifier()

Private Methods
---------------

.. method:: _lookup_loader(file_path: str) -> AbstractLoader

   Resolve the appropriate loader for a given file path. Classifies the file with
   ``RawScanClassifier.get_classification`` and looks the resulting
   ``RawScanType`` up in the registry.

   :param file_path: Path to the input file.
   :type file_path: str
   :returns: Loader instance corresponding to the file classification.
   :rtype: AbstractLoader
   :raises RuntimeError: If no loader is registered for the classification.

Public Methods
--------------

.. method:: load_file(file_path: str, loader: LoaderInterface | None = None) -> RawScan

   Load a single file into a :class:`RawScan`.

   If no loader is provided, one is inferred via classification.

   :param file_path: Path to the file.
   :type file_path: str
   :param loader: Optional explicit loader instance.
   :type loader: LoaderInterface, optional
   :returns: Loaded raw scan object.
   :rtype: RawScan

.. method:: load_files(file_paths: list[str], loader: LoaderInterface | None = None, quick: bool = Config["library.storage.raw.classification.quick"]) -> list[RawScan]

   Load multiple files into a list of :class:`RawScan` objects.

   When ``quick`` is enabled and no loader is provided, the loader is resolved
   once using the first file and reused for all subsequent files. This assumes
   homogeneous file types and improves performance.

   :param file_paths: List of file paths.
   :type file_paths: list[str]
   :param loader: Optional explicit loader instance.
   :type loader: LoaderInterface, optional
   :param quick: Enable single classification for batch loading.
   :type quick: bool
   :returns: List of loaded raw scans.
   :rtype: list[RawScan]

.. method:: load_folder(folder_path: str, loader: LoaderInterface | None = None, quick: bool = Config["library.storage.raw.classification.quick"]) -> list[RawScan]

   Load all files from a folder into :class:`RawScan` objects.

   Internally retrieves file paths from the filestore and delegates to
   :meth:`load_files`.

   :param folder_path: Path to the folder.
   :type folder_path: str
   :param loader: Optional explicit loader instance.
   :type loader: LoaderInterface, optional
   :param quick: Enable single classification for batch loading.
   :type quick: bool
   :returns: List of loaded raw scans.
   :rtype: list[RawScan]

Design Notes
------------

- **Singleton pattern** ensures a single shared controller instance.
- **Pluggable loaders** via ``LoaderRegistry`` support extensibility.
- **Performance optimization** via ``quick`` mode reduces redundant classification.
- Assumes consistent file types when ``quick=True``; misuse may lead to incorrect loader selection.
