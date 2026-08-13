Logging
=======

Overview
--------

TAVI uses Python's standard ``logging`` module configured via JSON, through
``tavi.meta.logging.init_logging()``.

Where init_logging() is called
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is **not** called from ``tavi/__init__.py`` and **not** called by the
``tavi`` console-script entry point (``tavi.__main__:execute``). The two call
sites are:

- ``tavi/__main__.py``, but only under ``if __name__ == "__main__"`` — i.e. when
  the app is started with ``python -m tavi`` (which is what ``pixi run tavi``
  does).
- ``tests/conftest.py``, so the test suite always has it configured.

Launching via the installed ``tavi`` script therefore leaves logging
unconfigured, and records fall back to Python's default handler-of-last-resort
(WARNING and above, no hostname, to stderr). Call ``init_logging()`` explicitly
if you need configured logging from an embedding process or a notebook:

.. code-block:: python

    from tavi.meta.logging import init_logging

    init_logging()

Configuration
-------------

The logging configuration is read with ``neutrons_standard.config.Resource`` from
``src/tavi/resources/logging_config.json`` and includes:

- **Formatters**: ``standard`` (hostname, timestamp, level, logger name, message)
  and ``simple`` (level and message only; defined but not currently attached to
  a handler)
- **Handlers**: a single ``console`` handler streaming to ``sys.stdout`` at INFO
- **Filters**: ``hostname_filter`` is injected into the config dict at runtime by
  ``init_logging`` rather than being declared in the JSON, because it needs the
  ``HostnameFilter`` class object
- **Root Logger**: ``""`` configured at INFO with ``propagate: true``

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

- Initialization: see `Where init_logging() is called`_ — it is not automatic for
  every entry point
- Hostname Filter: ``HostnameFilter`` is defined *inside* ``init_logging()`` and
  adds a ``hostname`` attribute (the short hostname, everything before the first
  dot) to each record. The ``standard`` formatter's ``%(hostname)s`` depends on
  it, so that formatter only works on handlers carrying the filter
- Configuration Source: ``src/tavi/resources/logging_config.json``, loaded via
  ``Resource.read`` and applied with ``logging.config.dictConfig``
