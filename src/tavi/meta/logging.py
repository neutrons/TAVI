"""Logging configuration for TAVI."""

import json
import logging.config
import socket

from neutrons_standard.config import Resource


def init_logging() -> None:
    """Initialize logging configuration with hostname filter."""

    class HostnameFilter(logging.Filter):
        """Add hostname to log records."""

        def filter(self, record: logging.LogRecord) -> bool:
            """Add hostname attribute to log record."""
            record.hostname = socket.gethostname().split(".")[0]
            return True

    LOGGING_CONFIG = json.loads(Resource.read("logging_config.json"))
    LOGGING_CONFIG["filters"] = {
        "hostname_filter": {
            # Reference the custom filter class
            "()": HostnameFilter,
        },
    }

    logging.config.dictConfig(LOGGING_CONFIG)
