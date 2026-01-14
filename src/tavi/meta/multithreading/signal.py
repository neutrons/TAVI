"""Signal."""

import asyncio
import inspect
import logging
import threading
import weakref
from typing import Any, Callable


class Signal:
    """PyQt-like implementation of cross thread callback communication."""

    def __init__(self, loop: asyncio.AbstractEventLoop) -> None:
        """Initialize signal with relevant event loop."""
        self._loop = loop
        self._loop_thread_id = None
        self._slots = []

    def bind_loop_thread(self) -> None:
        """Must be called from inside the running event loop thread."""
        self._loop_thread_id = threading.get_ident()

    def connect(self, slot: Callable) -> None:
        """Register a consumer to be emitted to later."""
        if inspect.ismethod(slot):
            self._slots.append(weakref.WeakMethod(slot))
        else:
            self._slots.append(weakref.ref(slot))

    def emit(self, *args: Any, **kwargs: Any) -> None:
        """Send data to Signal consumers."""
        # Determine connection type (Qt::AutoConnection semantics)
        is_loop_thread = threading.get_ident() == self._loop_thread_id

        for ref in list(self._slots):
            slot = ref()
            if slot is None:
                self._slots.remove(ref)
                continue

            if is_loop_thread:
                # Direct connection
                self._invoke_direct(slot, *args, **kwargs)
            else:
                # Queued connection
                self._invoke_queued(slot, *args, **kwargs)

    # --- internals ---

    def _invoke_direct(self, slot: Callable, *args: Any, **kwargs: Any) -> None:
        try:
            if inspect.iscoroutinefunction(slot):
                asyncio.create_task(slot(*args, **kwargs))
            else:
                slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (direct)")

    def _invoke_queued(self, slot: Callable, *args: Any, **kwargs: Any) -> None:
        if inspect.iscoroutinefunction(slot):
            self._loop.call_soon_threadsafe(
                asyncio.create_task,
                slot(*args, **kwargs),
            )
        else:
            self._loop.call_soon_threadsafe(
                self._safe_call,
                slot,
                *args,
                **kwargs,
            )

    @staticmethod
    def _safe_call(slot: Callable, *args: Any, **kwargs: Any) -> None:
        try:
            slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (queued)")
