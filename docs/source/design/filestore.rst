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

    def write_text_file(file_path: str, value: str) -> None
        """Write text content to a file."""

    def read_text_file(file_path: str) -> str
        """Read text content from a file."""

    def get_file_ext(file_path: str) -> str
        """Get the file extension."""

    def get_file_size_mb(file_path: str) -> float
        """Get file size in megabytes."""

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

- Writes to user application home directory
- Resolves tilde paths
- **Note**: Currently checks if value (content) starts with "/", not the path

Example:

.. code-block:: python

    filestore = LocalFileStore()
    # Writes to $USER_HOME/mydata/config.txt
    filestore.write_user_data_file("mydata/config.txt", "settings")

File Reading
~~~~~~~~~~~~

**read_text_file(file_path)**

- Reads entire file content as string
- Raises ``RuntimeError`` if file doesn't exist or is a directory
- Validates file existence before reading

Example:

.. code-block:: python

    content = filestore.read_text_file("/path/to/file.txt")

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

- Returns file extension including the dot (e.g., ".txt", ".json")
- Raises ``RuntimeError`` if file doesn't exist
- Works with multi-part extensions (e.g., ".tar.gz" returns ".gz")

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

1. Validates directory path exists
2. Raises ``RuntimeError`` if path is a file
3. Iterates through directory contents
4. Applies ``validate_file()`` to each item
5. Returns list of valid file paths (absolute)

Example:

.. code-block:: python

    valid_files = filestore.fetch_files_at("/data/scans")
    for filepath in valid_files:
        print(f"Processing: {filepath}")

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

Default value should be set in application configuration.

User Home Directory
~~~~~~~~~~~~~~~~~~~

User data files are written relative to:

.. code-block:: python

    Config["user.application.home"]  # typically ~/.<appname>

This supports tilde expansion for cross-platform compatibility.

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
        print(f"File error: {e}")  # "File at path ... does not exist"

    try:
        filestore.fetch_files_at("/data/file.txt")  # passing a file
    except RuntimeError as e:
        print(f"Path error: {e}")  # "Path ... is a file. A folder path is expected"

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

Use ``Config_override`` context manager for tests:

.. code-block:: python

    from tests.util.Config_helpers import Config_override

    with Config_override("library.filestore.raw.size-limit", 100):
        # Tests with custom size limit
        result = filestore.validate_file(filepath)

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
