"""Enumeration of the comparisons a filter criterion can make."""

from enum import StrEnum


class FilterOperator(StrEnum):
    """Comparison a ``FilterCriterion`` makes between a scan's metadata value and its own value."""

    IS = "is"
    IS_NOT = "is not"
    HAS = "has"
    HAS_NOT = "has not"
    EQUAL = "="
    NOT_EQUAL = "!="
    GREATER = ">"
    LESS = "<"
    GREATER_EQUAL = ">="
    LESS_EQUAL = "<="

    @property
    def is_numeric(self) -> bool:
        """Return whether this operator compares values as numbers rather than as text."""
        return self not in (FilterOperator.IS, FilterOperator.IS_NOT, FilterOperator.HAS, FilterOperator.HAS_NOT)
