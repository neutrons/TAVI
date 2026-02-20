"""Class each tavi base exception must inherit."""

import copy


class TaviError(Exception):
    """Base exception for tavi."""

    def __init__(self, message: str, stack_trace: str) -> None:
        """Initialize TaviError."""
        super().__init__(message)
        self.message = message
        self.stack_trace = stack_trace

    def __deepcopy__(self, memo: dict) -> "TaviError":
        """Create a new instance with the required arguments."""
        new_instance = type(self)(copy.deepcopy(self.message, memo), copy.deepcopy(self.stack_trace, memo))
        # Store the new instance in memo to handle circular references
        memo[id(self)] = new_instance
        return new_instance
