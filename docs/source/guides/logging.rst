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

Implementation Details
----------------------

- Initialization: ``init_logging()`` is called automatically in ``tavi/__init__.py``
- Hostname Filter: Custom filter adds ``hostname`` attribute to each log record
- Configuration Source: ``tavi.meta.logging.HostnameFilter`` and ``logging_config.json``
