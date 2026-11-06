"""Enumeration of the data scopes a filter criterion can select."""

from enum import StrEnum


class FilterScope(StrEnum):
    """Key set and grouping selected by a criterion for matching."""

    META = "metadata"
    COLUMN_FULL = "col full"
    COLUMN_PARTIAL = "col part"
