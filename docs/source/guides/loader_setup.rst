Setting Up a New Loader
=======================

The TAVI framework provides a plugin-based architecture for handling different scan file formats through the **Loader** system. This guide explains how to create and register a new loader so it gets automatically picked up by the `LoaderRegistry` and integrated into the `RawScanClassifier` flow.

Overview
--------

The loader system consists of three main components:

1. **LoaderInterface**: Abstract interface that all loaders must implement
2. **LoaderRegistry**: Singleton that manages and coordinates all loaders
3. **RawScanClassifier**: Uses the registry to automatically classify scan file types

The classification flow works as follows:

.. code-block:: text

    File Path → RawScanClassifier → LoaderRegistry → All Loaders
                                         ↓
                        Each loader scores the file
                                         ↓
                        Highest scoring loader wins
                                         ↓
                                   File Type Identified

Step 1: Create Your Loader Class
--------------------------------

Create a new Python file in ``src/tavi/library/storage/loader/`` directory. Your loader must inherit from ``AbstractLoader``.

Example: Creating a loader for a hypothetical "MyFormat" scan type

.. code-block:: python

    """MyFormat scan file loader."""

    from typing import Any

    from tavi.library.data.enum.raw_scan_type import RawScanType
    from tavi.library.data.scan import Scan, ScanData, ScanMetadata
    from tavi.library.storage.interface.file_store_interface import FileStoreInterface
    from tavi.library.storage.loader.interface.base import AbstractLoader


    class MyFormatLoader(AbstractLoader):
        """Loader for MyFormat scan files."""

        def __init__(self, filestore: FileStoreInterface) -> None:
            """Initialize MyFormat loader."""
            super().__init__(filestore)

        def load(self, path: str) -> Scan:
            """Load scan data from file."""
            # Implement your file loading logic
            pass

        def get_scan_type(self) -> RawScanType:
            """Get scan type identifier."""
            # Must match the enum value you add in Step 2
            return RawScanType.MyFormat

        def get_score(self, path: str) -> float:
            """
            Score how confident this loader is for a given file.

            Return a score between 0.0 and 1.0 where:
            - 0.0: This loader cannot handle this file
            - 0.1-0.5: Maybe this loader can handle it
            - 0.51-1.0: Strongly confident this loader can handle it

            The loader with the highest score will be selected.
            """
            # Implement your file format detection logic
            if path.endswith('.myformat'):
                return 1.0  # Strong match
            return 0.0

        def parse_metadata(self, path: str) -> ScanMetadata:
            """Parse metadata from the file."""
            pass

        def parse_scan_values(self, path: str) -> ScanData:
            """Parse scan data values from the file."""
            pass

        def parse_external_metadata(self, path: str) -> dict[str, Any]:
            """Parse any external metadata associated with the file."""
            pass

        def adapt_scan_data(self, meta: ScanMetadata, values: ScanData) -> Scan:
            """Combine metadata and values into a Scan object."""
            pass

Step 2: Add Enum Value for Your Format
---------------------------------------

Add your format to the ``RawScanType`` enum in ``src/tavi/library/data/enum/raw_scan_type.py``:

.. code-block:: python

    """Enumeration of raw scan types."""

    from enum import StrEnum


    class RawScanType(StrEnum):
        """Enumeration of supported raw scan file types."""

        ORNLSpice = "ORNLSpice"
        MyFormat = "MyFormat"  # Add your format here
        NONE = "None"

Step 3: Register Your Loader
-----------------------------

Register your loader in the ``LoaderRegistry`` singleton located at ``src/tavi/library/storage/loader/loader_registry.py``:

.. code-block:: python

    """Loader registry singleton."""

    from neutrons_standard.decorators.singleton import Singleton

    from tavi.library.storage.interface.filestore_interface import Filestore
    from tavi.library.storage.loader.default_loader import DefaultLoader
    from tavi.library.storage.loader.interface.base import AbstractLoader
    from tavi.library.storage.loader.my_format_loader import MyFormatLoader  # Add import
    from tavi.library.storage.loader.ornl_spice_loader import ORNLSpiceLoader


    @Singleton
    class LoaderRegistry:
        """Registry for managing loaders."""

        def __init__(self, filestore: Filestore) -> None:
            """Initialize registry with filestore."""
            self.registry: dict[str, AbstractLoader] = {}
            self.set_filestore(filestore)

            # Register loaders in order of priority (highest first)
            self._register_loader(ORNLSpiceLoader(self.filestore))
            self._register_loader(MyFormatLoader(self.filestore))  # Add your loader
            self._register_loader(DefaultLoader(self.filestore))

        # ... rest of the implementation

**Important**: Register more specific loaders **before** more general ones. The `DefaultLoader` should always be last as it returns a score of 0 for everything.

Step 4: Test Your Loader
------------------------

Create tests to verify your loader:

1. **Scoring**: Test that `get_score()` returns appropriate scores for your format
2. **Loading**: Test that `load()` correctly reads and parses files
3. **Classification**: Test that `RawScanClassifier` correctly identifies your file type

Example test:

.. code-block:: python

    from tavi.backend.classification.raw_scan_classifier import RawScanClassifier
    from tavi.library.data.enum.raw_scan_type import RawScanType


    def test_myformat_classification():
        classifier = RawScanClassifier()
        result = classifier.get_classification("path/to/file.myformat")
        assert result == RawScanType.MyFormat


How Classification Works
------------------------

When ``RawScanClassifier.get_classification(file_path)`` is called:

1. It retrieves all registered loaders from `LoaderRegistry`
2. It calls `get_score(file_path)` on **each** loader
3. It tracks which loader returned the highest score
4. It returns the `get_scan_type()` of the winning loader
5. If no loader scores above 0, it returns `RawScanType.NONE` (from DefaultLoader)

.. code-block:: python

    # Example with multiple loaders
    classifier = RawScanClassifier()

    # LoaderRegistry contains: [ORNLSpiceLoader, MyFormatLoader, DefaultLoader]
    # For "scan.myformat":
    #   ORNLSpiceLoader.get_score("scan.myformat") → 0.0
    #   MyFormatLoader.get_score("scan.myformat") → 1.0
    #   DefaultLoader.get_score("scan.myformat") → 0.0
    # Winner: MyFormatLoader
    file_type = classifier.get_classification("scan.myformat")
    # Result: RawScanType.MyFormat

Best Practices
--------------

1. **Implement Robust Scoring**: Your `get_score()` method should be fast and use multiple heuristics, returning a float between 0.0 and 1.0:

   - Check file extension
   - Validate magic bytes (file header)
   - Check for format-specific markers
   - Examine file structure

2. **Handle Edge Cases**: Return 0.0 for files you can't handle, not negative scores

3. **Document Your Format**: Add docstrings explaining what file formats your loader supports

4. **Fail Gracefully**: If file parsing fails in `load()`, raise clear exceptions with context

5. **Keep Score Ranges Consistent**:

   - 0.0 = definitely not your format
   - 0.01-0.5 = uncertain
   - 0.51-1.0 = confident match

Example: Real-World Implementation
-----------------------------------

The ``ORNLSpiceLoader`` demonstrates a production-ready implementation:

.. code-block:: python

    class ORNLSpiceLoader(AbstractLoader):
        """Loader for ORNL Spice format scan files."""

        def __init__(self, filestore: FileStoreInterface) -> None:
            """Initialize ORNL Spice loader with classifier."""
            super().__init__(filestore)
            self.classifier = RuleBasedClassifier()
            self.classification_rules = ORNLSpiceRuleSet()

        def get_score(self, path: str) -> int:
            """Get score for scan using rule-based classification."""
            # Uses a dedicated RuleBasedClassifier for intelligent scoring
            return self.classifier.get_score(path, self.classification_rules)

This example shows how you can:
- Use helper classifiers for complex format detection
- Leverage rule sets for sophisticated scoring logic
- Maintain clean separation of concerns

Troubleshooting
---------------

**My loader never gets selected**
    - Check that `get_score()` returns a higher score than other loaders for your file type (between 0.0 and 1.0)
    - Verify your loader is registered in `LoaderRegistry.__init__()`
    - Ensure `get_scan_type()` returns the correct enum value

**RawScanType.NONE is returned**
    - All loaders returned score 0.0, or no loaders were registered
    - Verify your loader's `get_score()` implementation returns values between 0.0 and 1.0
    - Check that your `RawScanType` enum value exists

**Precommit checks fail**
    - Add module and class docstrings to your loader file
    - Add docstrings to all public methods
    - Run `ruff format` to auto-fix formatting issues

See Also
--------

- :doc:`rule_based_classifier` - For implementing complex classification logic
- Loader Interface: ``src/tavi/library/storage/loader/interface/loader_interface.py``
- Existing Implementations: ``src/tavi/library/storage/loader/ornl_spice_loader.py``
