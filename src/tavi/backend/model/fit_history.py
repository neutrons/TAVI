"""Per-member undo/redo history of a fit's saved members."""

from collections import deque
from collections.abc import Callable, Container
from typing import Optional

from tavi.library.data.fit_entry import FitMember
from tavi.library.data.scan import UUID

# Undo is for dialing one member in, not an audit trail - a handful of steps back is plenty, and
# anything older is dropped first.
FIT_HISTORY_LIMIT = 10

FitMemberKey = tuple[UUID, UUID]
"""(fit uuid, source scan uuid) - a member's identity inside its fit is its source scan."""


class FitHistory:
    """Undo/redo stacks of earlier saved states, one pair per member of a fit."""

    def __init__(self, limit: int = FIT_HISTORY_LIMIT) -> None:
        """Init empty stacks, each keeping at most ``limit`` states."""
        self._limit = limit
        self._undo: dict[FitMemberKey, deque[FitMember]] = {}
        self._redo: dict[FitMemberKey, deque[FitMember]] = {}

    def record(self, fit_uuid: UUID, previous: FitMember) -> None:
        """Remember the state a fresh save is about to replace; a new step forks history, so redo is cleared."""
        key = (fit_uuid, previous.source_scan_uuid)
        self._stack(self._undo, key).append(previous)
        self._redo.pop(key, None)

    def undo(self, fit_uuid: UUID, current: FitMember) -> Optional[FitMember]:
        """Pop the state before ``current``, keeping ``current`` for redo - ``None`` if there's nothing to undo."""
        return self._step(fit_uuid, current, self._undo, self._redo)

    def redo(self, fit_uuid: UUID, current: FitMember) -> Optional[FitMember]:
        """Pop the state after ``current``, keeping ``current`` for undo - ``None`` if there's nothing to redo."""
        return self._step(fit_uuid, current, self._redo, self._undo)

    def can_undo(self, key: FitMemberKey) -> bool:
        """Whether ``key``'s member has an earlier state to go back to."""
        return bool(self._undo.get(key))

    def can_redo(self, key: FitMemberKey) -> bool:
        """Whether ``key``'s member has an undone state to go forward to."""
        return bool(self._redo.get(key))

    def forget(self, fit_uuids: Container[UUID] = (), scan_uuids: Container[UUID] = ()) -> None:
        """Drop the history of every member of ``fit_uuids``, and of every member fit against ``scan_uuids``."""
        self._drop(lambda key: key[0] in fit_uuids or key[1] in scan_uuids)

    def _step(
        self,
        fit_uuid: UUID,
        current: FitMember,
        source: dict[FitMemberKey, deque[FitMember]],
        target: dict[FitMemberKey, deque[FitMember]],
    ) -> Optional[FitMember]:
        key = (fit_uuid, current.source_scan_uuid)
        stack = source.get(key)
        if not stack:
            return None
        self._stack(target, key).append(current)
        return stack.pop()

    def _stack(self, stacks: dict[FitMemberKey, deque[FitMember]], key: FitMemberKey) -> deque[FitMember]:
        # maxlen makes the deque push the oldest state out once the limit is reached.
        return stacks.setdefault(key, deque(maxlen=self._limit))

    def _drop(self, matches: Callable[[FitMemberKey], bool]) -> None:
        for stacks in (self._undo, self._redo):
            for key in [key for key in stacks if matches(key)]:
                del stacks[key]
