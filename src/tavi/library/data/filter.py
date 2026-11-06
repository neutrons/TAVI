"""Scan filter domain models."""

import operator
from typing import Any, Callable

from neutrons_standard.config import Config
from pydantic import BaseModel, ConfigDict, Field

from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.enum.filter_scope import FilterScope
from tavi.library.data.scan import UUID, Scan, UUIDFactory

_NUMERIC_COMPARISONS: dict[FilterOperator, Callable[[float, float], bool]] = {
    FilterOperator.EQUAL: operator.eq,
    FilterOperator.NOT_EQUAL: operator.ne,
    FilterOperator.GREATER: operator.gt,
    FilterOperator.LESS: operator.lt,
    FilterOperator.GREATER_EQUAL: operator.ge,
    FilterOperator.LESS_EQUAL: operator.le,
}


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_str(value: Any) -> str:
    return str(value).strip().casefold()


class FilterCriterion(BaseModel):
    """One atomic statement of a filter, comparing a scan's metadata entry ``key`` against ``value``."""

    scope: FilterScope
    key: str
    operator: FilterOperator
    value: str

    def validate_criterion(self) -> list[str]:
        """Return why this criterion cannot be evaluated, or an empty list if it can."""
        errors = []
        if not self.key.strip():
            errors.append("Criterion has no key.")
        if not self.value.strip():
            errors.append(f"Criterion on {self.key!r} has no value.")
        elif self.operator.is_numeric and _as_float(self.value) is None:
            errors.append(
                f"Criterion on {self.key!r} compares with {self.operator!s} but {self.value!r} is not a number."
            )
        if not self.operator.is_numeric and self.scope in (FilterScope.COLUMN_FULL, FilterScope.COLUMN_PARTIAL):
            errors.append(f"Columnar criterion compares with {self.operator!s} but only supports numeric comparison.")
        return errors

    def matches(self, scan: Scan) -> bool:
        """
        Return whether ``scan`` satisfies this criterion.

        A scan without ``key`` never matches, whatever the operator - "is not" should not sweep in
        every scan from an instrument that doesn't record the field at all. Text comparisons ignore
        case; a numeric comparison against a value that isn't a number does not match.
        """
        match self.scope:
            case FilterScope.META:
                return self._match_metadata(scan.metadata.data)
            case FilterScope.COLUMN_FULL:
                return self._match_columnar(scan.data.data, partial=False)
            case FilterScope.COLUMN_PARTIAL:
                return self._match_columnar(scan.data.data, partial=True)
        return False

    def _match_metadata(self, metadata: dict[str, Any]) -> bool:
        """Return whether ``metadata`` satisfies this criterion."""
        if self.key not in metadata:
            return False
        meta_value = metadata.get(self.key, None)

        if self.operator.is_numeric:
            actual_number, expected_number = _as_float(meta_value), _as_float(self.value)
            if actual_number is None or expected_number is None:
                return False
            return _NUMERIC_COMPARISONS[self.operator](actual_number, expected_number)

        actual_text, expected_text = _as_str(meta_value), _as_str(self.value)
        match self.operator:
            case FilterOperator.IS:
                return actual_text == expected_text
            case FilterOperator.IS_NOT:
                return actual_text != expected_text
            case FilterOperator.HAS:
                return expected_text in actual_text
            case FilterOperator.HAS_NOT:
                return expected_text not in actual_text
        return False

    def _match_columnar(self, columnar_data: dict[str, list[Any]], partial: bool = False) -> bool:
        """Return whether ``column data`` satisfies this criterion."""
        if self.key not in columnar_data:
            return False
        col_data = columnar_data.get(self.key, [None])
        _logical_join: Callable = any if partial else all

        if self.operator.is_numeric:
            col_numbers, expected_number = list(map(_as_float, col_data)), _as_float(self.value)
            if expected_number is None or _logical_join(num is None for num in col_numbers):
                return False
            return _logical_join(
                _NUMERIC_COMPARISONS[self.operator](num, expected_number) for num in col_numbers if num is not None
            )

        # Made assuming columnar data can also support strings
        # Left in case columnar data supports strings in the future
        # col_text, expected_text = list(map(_as_str, col_data)), _as_str(self.value)
        # match self.operator:
        #     case FilterOperator.IS:
        #         return _logical_join(text == expected_text for text in col_text)
        #     case FilterOperator.IS_NOT:
        #         return _logical_join(text != expected_text for text in col_text)
        #     case FilterOperator.HAS:
        #         return _logical_join(expected_text in text for text in col_text)
        #     case FilterOperator.HAS_NOT:
        #         return _logical_join(expected_text not in text for text in col_text)
        return False


class ScanFilter(BaseModel):
    """A named, saved set of criteria; a scan passes when it satisfies every one of them."""

    uuid: UUID = UUIDFactory()
    name: str
    criteria: list[FilterCriterion] = Field(default_factory=list)
    criteria_union: bool = True

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def validate_filter(self) -> list[str]:
        """Return why this filter cannot be applied, or an empty list if it can."""
        errors = [] if self.name.strip() else ["Filter has no name."]
        for criterion in self.criteria:
            errors.extend(criterion.validate_criterion())
        return errors

    def matches(self, scan: Scan) -> bool:
        """Return whether ``scan``'s metadata satisfies every criterion; a filter without criteria passes all."""
        if self.criteria_union:
            return all(criterion.matches(scan) for criterion in self.criteria)
        return any(criterion.matches(scan) for criterion in self.criteria)


def sample_filter_uuid() -> UUID:
    """Return the uuid of the built-in Sample filter, as set in ``filter.sample.uuid``."""
    return UUID(value=Config["filter.sample.uuid"])


def sample_filter() -> ScanFilter:
    """Return the built-in, read-only Sample filter every project starts with, as set under ``filter.sample``."""
    criterion = FilterCriterion(
        scope=FilterScope(Config["filter.sample.criterion.scope"]),
        key=Config["filter.sample.criterion.key"],
        operator=FilterOperator(Config["filter.sample.criterion.operator"]),
        value=str(Config["filter.sample.criterion.value"]),
    )
    return ScanFilter(uuid=sample_filter_uuid(), name=Config["filter.sample.name"], criteria=[criterion])
