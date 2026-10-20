"""Scan filter domain models."""

import operator
from typing import Any, Callable

from neutrons_standard.config import Config
from pydantic import BaseModel, ConfigDict, Field

from tavi.library.data.enum.filter_operator import FilterOperator
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


class FilterCriterion(BaseModel):
    """One atomic statement of a filter, comparing a scan's metadata entry ``key`` against ``value``."""

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
        return errors

    def matches(self, metadata: dict[str, Any]) -> bool:
        """
        Return whether ``metadata`` satisfies this criterion.

        A scan without ``key`` never matches, whatever the operator - "is not" should not sweep in
        every scan from an instrument that doesn't record the field at all. Text comparisons ignore
        case; a numeric comparison against a value that isn't a number does not match.
        """
        if self.key not in metadata:
            return False
        actual = metadata[self.key]

        if self.operator.is_numeric:
            actual_number, expected_number = _as_float(actual), _as_float(self.value)
            if actual_number is None or expected_number is None:
                return False
            return _NUMERIC_COMPARISONS[self.operator](actual_number, expected_number)

        actual_text, expected_text = str(actual).strip().casefold(), self.value.strip().casefold()
        match self.operator:
            case FilterOperator.IS:
                return actual_text == expected_text
            case FilterOperator.IS_NOT:
                return actual_text != expected_text
            case FilterOperator.HAS:
                return expected_text in actual_text
            case FilterOperator.HAS_NOT:
                return expected_text not in actual_text


class ScanFilter(BaseModel):
    """A named, saved set of criteria; a scan passes when it satisfies every one of them."""

    uuid: UUID = UUIDFactory()
    name: str
    criteria: list[FilterCriterion] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def validate_filter(self) -> list[str]:
        """Return why this filter cannot be applied, or an empty list if it can."""
        errors = [] if self.name.strip() else ["Filter has no name."]
        for criterion in self.criteria:
            errors.extend(criterion.validate_criterion())
        return errors

    def matches(self, scan: Scan) -> bool:
        """Return whether ``scan``'s metadata satisfies every criterion; a filter without criteria passes all."""
        return all(criterion.matches(scan.metadata.data) for criterion in self.criteria)


def sample_filter_uuid() -> UUID:
    """Return the uuid of the built-in Sample filter, as set in ``filter.sample.uuid``."""
    return UUID(value=Config["filter.sample.uuid"])


def sample_filter() -> ScanFilter:
    """Return the built-in, read-only Sample filter every project starts with, as set under ``filter.sample``."""
    criterion = FilterCriterion(
        key=Config["filter.sample.criterion.key"],
        operator=FilterOperator(Config["filter.sample.criterion.operator"]),
        value=str(Config["filter.sample.criterion.value"]),
    )
    return ScanFilter(uuid=sample_filter_uuid(), name=Config["filter.sample.name"], criteria=[criterion])
