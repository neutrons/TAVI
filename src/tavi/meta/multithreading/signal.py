"""Signal."""

import asyncio
import inspect
import logging
import threading
import weakref
from typing import Any, Callable


class Signal:
    """PyQt-like implementation of cross-thread callback communication."""

    def __init__(self, loop: asyncio.AbstractEventLoop) -> None:
        """Initialize signal with relevant event loop."""
        self._loop = loop
        self._loop_thread_id: int | None = None
        self._slots: list[weakref.ReferenceType] = []

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
        if self._loop_thread_id is None:
            raise RuntimeError("Signal.bind_loop_thread() was never called.")

        is_loop_thread = threading.get_ident() == self._loop_thread_id

        for ref in list(self._slots):
            slot = ref()
            if slot is None:
                self._slots.remove(ref)
                continue

            if inspect.iscoroutinefunction(slot):
                # ALWAYS schedule coroutine slots
                self._schedule_coroutine(slot, *args, **kwargs)
            else:
                if is_loop_thread:
                    self._invoke_direct(slot, *args, **kwargs)
                else:
                    self._invoke_queued(slot, *args, **kwargs)

    # --- internals ---

    def _invoke_direct(self, slot: Callable, *args: Any, **kwargs: Any) -> None:
        try:
            slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (direct)")

    def _invoke_queued(self, slot: Callable, *args: Any, **kwargs: Any) -> None:
        self._loop.call_soon_threadsafe(
            self._safe_call,
            slot,
            *args,
            **kwargs,
        )

    def _schedule_coroutine(self, slot: Callable, *args: Any, **kwargs: Any) -> None:
        coro = slot(*args, **kwargs)

        def _create_task() -> None:
            try:
                asyncio.create_task(coro)
            except Exception:
                logging.exception("Signal coroutine slot failed")

        self._loop.call_soon_threadsafe(_create_task)

    @staticmethod
    def _safe_call(slot: Callable, *args: Any, **kwargs: Any) -> None:
        try:
            slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (queued)")
