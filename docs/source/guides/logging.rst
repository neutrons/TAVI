Logging
=======

Overview
--------

TAVI uses Python's standard ``logging`` module configured via JSON. Logging is automatically initialized when the application starts via ``tavi.meta.logging.init_logging()``.

Configuration
-------------

The logging configuration is defined in ``tavi/resources/logging_config.json`` and includes:

- **Formatters**: Define output format with hostname, timestamp, level, logger name, and message
- **Handlers**: Console handler streams to stdout at INFO level
- **Filters**: Custom hostname filter adds system hostname to all log records
- **Root Logger**: Configured at INFO level with propagation enabled

Log Format
~~~~~~~~~~

.. code-block:: text

    hostname - YYYY-MM-DD HH:MM:SS - LEVEL    - logger.name - message

Example:

.. code-block:: text

    my-server - 2026-03-09 13:10:54 - INFO     - tavi.library.storage - Loading scan data

Using Logging in Code
---------------------

Get a logger for your module:

.. code-block:: python

    import logging

    logger = logging.getLogger(__name__)

    logger.info("Processing scan data")
    logger.warning("Unusual value detected")
    logger.error("Failed to load file")

Changing Log Levels
~~~~~~~~~~~~~~~~~~~

To adjust logging verbosity, modify the ``level`` in ``logging_config.json``:

- ``DEBUG``: Detailed diagnostic information
- ``INFO``: General informational messages (default)
- ``WARNING``: Warning messages
- ``ERROR``: Error messages only
- ``CRITICAL``: Critical errors only

Logging Guidelines
------------------

Choose the appropriate log level based on context:

**DEBUG**
    Internal implementation details, variable values, function entry/exit. Use during development and troubleshooting.

    .. code-block:: python

        logger.debug(f"Loading file: {filepath}")
        logger.debug(f"Parsed {count} data points")

**INFO**
    General application flow and status. User-facing information that indicates normal operation.

    .. code-block:: python

        logger.info("Scan data loaded successfully")
        logger.info("Processing completed in 2.5 seconds")

**WARNING**
    Unexpected but recoverable situations. Data validation issues, deprecated usage, unusual conditions.

    .. code-block:: python

        logger.warning("Missing optional metadata field")
        logger.warning("Scan data appears incomplete")

**ERROR**
    Serious failures that prevent operations but don't crash the application.

    .. code-block:: python

        logger.error("Failed to load file: invalid format")
        logger.error("Database connection lost")

**CRITICAL**
    Severe errors that may require application shutdown or manual intervention.

    .. code-block:: python

        logger.critical("Configuration file missing or corrupted")
        logger.critical("Unable to initialize core service")

Implementation Details
----------------------

- Initialization: ``init_logging()`` is called automatically in ``tavi/__init__.py``
- Hostname Filter: Custom filter adds ``hostname`` attribute to each log record
- Configuration Source: ``tavi.meta.logging.HostnameFilter`` and ``logging_config.json``
