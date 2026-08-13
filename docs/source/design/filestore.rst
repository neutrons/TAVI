File Storage System
===================

The file storage system provides an abstraction for managing file operations,
with a local filesystem implementation. This design allows for flexible storage
backends while maintaining a consistent interface for file handling, validation,
and user data management.

Core Concepts
-------------

The storage system distinguishes between two primary concerns:

- **Interface** (``FileStoreInterface``)
  Defines the contract for any file storage implementation.

- **Implementation** (``LocalFileStore``)
  Concrete implementation for local filesystem operations.

Architecture Overview
---------------------

Components:

- ``FileStoreInterface`` - Abstract base class defining storage contract
- ``LocalFileStore`` - Concrete implementation for local filesystem
- ``Config`` - Configuration system for runtime settings
- ``Path`` (pathlib) - Cross-platform path handling

Responsibilities by method:

- **File Reading/Writing**: Read and write text files with automatic creation
- **File Validation**: Check file existence, type, and size constraints
- **File Inspection**: Extract metadata (extension, size)
- **Directory Scanning**: Fetch and filter files from directories
- **User Data Storage**: Write application data to user-specific locations

FileStoreInterface
------------------

Abstract base class defining the storage contract.

Methods:

.. code-block:: python

    def fetch_files_at(path: str) -> list[str]
        """Fetch all valid files from a directory."""

    def validate_file(file_path: str) -> bool
        """Validate if a file meets storage requirements."""

    def write_user_data_file(file_subpath: str, value: str) -> None
        """Write application data to user home directory."""

    def read_user_data_file(file_subpath: str) -> str
        """Read application data from user home directory."""

    def write_text_file(file_path: str, value: str) -> None
        """Write text content to a file."""

    def read_text_file(file_path: str) -> str
        """Read text content from a file."""

    def get_file_ext(file_path: str) -> str
        """Get the file extension."""

    def get_file_size_mb(file_path: str) -> float
        """Get file size in megabytes."""

    def get_file_name(file_path: str) -> str
        """Get the file name."""

    def get_parent(file_path: str) -> str
        """Get the parent directory of a path."""

    def join_path(root_path: str, target_path: str) -> str
        """Join two paths."""

LocalFileStore Implementation
-----------------------------

Concrete implementation for local filesystem operations.

File Writing
~~~~~~~~~~~~

**write_text_file(file_path, value)**

- Creates file if it doesn't exist
- Overwrites existing files
- Creates parent directories as needed (if using Path operations elsewhere)
- Logs all write operations

Example:

.. code-block:: python

    filestore = LocalFileStore()
    filestore.write_text_file("/path/to/file.txt", "content")

**write_user_data_file(file_subpath, value)**

- Writes under the user application home directory, ``Config["user.application.home"]``
  (``~/.TAVI/`` by default)
- Expands the tilde and creates the home directory if it is missing
- Strips a leading ``/`` from ``file_subpath``, so an absolute-looking subpath is
  still resolved relative to the user home rather than the filesystem root
- Delegates the actual write to ``write_text_file``

Example:

.. code-block:: python

    filestore = LocalFileStore()
    # Writes to ~/.TAVI/mydata/config.txt
    filestore.write_user_data_file("mydata/config.txt", "settings")

File Reading
~~~~~~~~~~~~

**read_text_file(file_path)**

- Reads entire file content as a UTF-8 string
- Raises ``RuntimeError`` if the path doesn't exist or is not a file
- Validates via ``_is_real_file(..., throws=True)`` before reading

Example:

.. code-block:: python

    content = filestore.read_text_file("/path/to/file.txt")

**read_user_data_file(file_subpath)**

- Mirror of ``write_user_data_file``: resolves ``file_subpath`` under the user
  application home, then delegates to ``read_text_file``
- Raises ``RuntimeError(f"File does not exist: {user_data_file}")`` when missing

Example:

.. code-block:: python

    raw_settings = filestore.read_user_data_file("settings.yaml")

File Validation & Inspection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**validate_file(file_path)**

Checks if a file meets storage requirements:

1. File must exist
2. Path must point to a file (not directory)
3. File size must not exceed ``library.filestore.raw.size-limit`` (MB)

Returns ``True`` if all checks pass, ``False`` otherwise.

Example:

.. code-block:: python

    if filestore.validate_file("/path/to/data.dat"):
        process_file("/path/to/data.dat")

**get_file_ext(file_path)**

- Returns file extension including the dot (e.g., ".txt", ".json"), via ``Path.suffix``
- Raises ``RuntimeError`` if file doesn't exist
- Works with multi-part extensions (e.g., ".tar.gz" returns ".gz")

**get_file_name(file_path)**

- Returns the final path component including the extension (``Path.name``)
- Raises ``RuntimeError`` if file doesn't exist
- Used by ``InstrumentInFilenameRule`` to test the leading token of a SPICE filename

**get_file_size_mb(file_path)**

- Returns file size in megabytes (float)
- Raises ``RuntimeError`` if file doesn't exist
- Conversion: ``size_bytes / (1024 * 1024)``

Example:

.. code-block:: python

    size_mb = filestore.get_file_size_mb("/path/to/data.h5")
    if size_mb > 1000:
        print(f"Large file: {size_mb:.2f} MB")

Directory Operations
~~~~~~~~~~~~~~~~~~~~

**fetch_files_at(path)**

Retrieves and filters files from a directory.

Process:

1. Raises ``RuntimeError`` if the path does not exist
2. Raises ``RuntimeError`` if path is a file
3. Iterates through directory contents (non-recursive, ``Path.iterdir``)
4. Applies ``validate_file()`` to each item
5. Returns a **sorted** list of valid absolute file paths

The sort is what makes folder loading deterministic — the order scans appear in
the project tree matches the lexical order of their filenames.

Example:

.. code-block:: python

    valid_files = filestore.fetch_files_at("/data/scans")
    for filepath in valid_files:
        print(f"Processing: {filepath}")

Path Helpers
~~~~~~~~~~~~

``get_parent(file_path)`` and ``join_path(root_path, target_path)`` are thin
``pathlib`` wrappers that do **not** touch the filesystem, so they work on paths
that do not exist yet. ``ORNLSpiceLoader`` uses them to resolve a scan's sibling
``UBConf`` directory from the ``.dat`` file's own path.

Internal Methods
~~~~~~~~~~~~~~~~

**_is_real_file(file_path, throws=False)**

Helper method to validate if path is a real file.

Parameters:

- ``file_path`` - Path object to check
- ``throws`` - If True, raises ``RuntimeError`` on failure

Returns:

- ``True`` if path exists and is a file
- ``False`` if path doesn't exist or is not a file
- Raises ``RuntimeError`` if ``throws=True`` and validation fails

Configuration
-------------

Size Limit
~~~~~~~~~~

The maximum file size for validation is controlled by:

.. code-block:: python

    Config["library.filestore.raw.size-limit"]  # in MB

The default lives in ``src/tavi/resources/application.yml`` and is currently
**1 MB**. Any file larger than that fails ``validate_file()`` and is silently
skipped by ``fetch_files_at``.

User Home Directory
~~~~~~~~~~~~~~~~~~~

User data files are written relative to:

.. code-block:: python

    Config["user.application.home"]  # ~/.TAVI/ by default

The value is expanded with ``Path.expanduser()``, so a leading ``~`` resolves per
platform.

Error Handling
--------------

Exception Types
~~~~~~~~~~~~~~~

``RuntimeError`` is raised for:

- File operations on nonexistent files
- Directory operations on file paths
- Directory operations on nonexistent paths
- Invalid file types during validation

Common Scenarios:

.. code-block:: python

    try:
        content = filestore.read_text_file("/missing/file.txt")
    except RuntimeError as e:
        # "File at path `/missing/file.txt` does not exist or is not a file."
        print(f"File error: {e}")

    try:
        filestore.fetch_files_at("/data/file.txt")  # passing a file
    except RuntimeError as e:
        # "Path `/data/file.txt` is a file.  A folder path is expected."
        print(f"Path error: {e}")

Validation Failures
~~~~~~~~~~~~~~~~~~~

``validate_file()`` returns ``False`` instead of raising exceptions:

- File doesn't exist
- Path is a directory
- File exceeds size limit

This allows graceful filtering without exception handling:

.. code-block:: python

    # Graceful filtering
    files = filestore.fetch_files_at("/data")
    # Returns only files that pass all validation checks

Usage Examples
--------------

Reading Scan Data
~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.storage.local_file_store import LocalFileStore

    filestore = LocalFileStore()

    # Fetch all valid scan files
    scan_files = filestore.fetch_files_at("/experiments/exp123/scans")

    for filepath in scan_files:
        # Check before processing
        if filestore.validate_file(filepath):
            ext = filestore.get_file_ext(filepath)
            size = filestore.get_file_size_mb(filepath)
            content = filestore.read_text_file(filepath)
            process_scan(content, size)

Storing User Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    filestore = LocalFileStore()

    config_data = """
    [settings]
    theme=dark
    resolution=1920x1080
    """

    # Automatically writes to user home
    filestore.write_user_data_file(
        "config/settings.ini",
        config_data
    )

Testing Considerations
----------------------

Test Fixtures
~~~~~~~~~~~~~

Use ``TemporaryDirectory`` for isolated tests:

.. code-block:: python

    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as tmpdir:
        filestore.write_text_file(f"{tmpdir}/test.txt", "data")
        # Directory auto-cleaned after context

Configuration Overrides
~~~~~~~~~~~~~~~~~~~~~~~

Use the ``Config_override`` context manager from ``tests/util/Config_helpers.py``:

.. code-block:: python

    from util.Config_helpers import Config_override

    with Config_override("library.filestore.raw.size-limit", 100):
        # Tests with custom size limit
        result = filestore.validate_file(filepath)

``tests/conftest.py`` also exposes it as the ``Config_override_fixture`` fixture,
which unwinds every override at teardown:

.. code-block:: python

    def test_large_file(Config_override_fixture):
        Config_override_fixture("library.filestore.raw.size-limit", 100)
        assert filestore.validate_file(filepath)

Design Characteristics
----------------------

- **Interface-based**: Abstraction allows multiple implementations
- **Error-aware**: Clear separation between validation failure and exceptions
- **Configuration-driven**: Runtime settings without code changes
- **Path-safe**: Uses pathlib for cross-platform compatibility
- **Simple semantics**: Direct mapping to filesystem operations
- **Logging**: All writes are logged for audit trails

Recommended Practices
---------------------

- Validate files before processing: Use ``validate_file()`` to filter
- Handle size constraints: Check ``get_file_size_mb()`` for large files
- Use user data API: Call ``write_user_data_file()`` for application config
- Log operations: Leverage built-in logging for debugging
- Test with temp directories: Ensure tests clean up after themselves
- Override config in tests: Use ``Config_override`` for testability
- Catch RuntimeError: Be specific about which operations may fail

Future Considerations
---------------------

- Remote storage backends (S3, cloud storage)
- Async file operations
- File compression support
- Additional metadata extraction (creation time, permissions)
- Streaming for large files
- File caching layer
