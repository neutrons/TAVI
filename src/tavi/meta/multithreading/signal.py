import asyncio
import threading
import inspect
import weakref
import logging


class Signal:
    def __init__(self, loop: asyncio.AbstractEventLoop, *args):
        self._loop = loop
        self._loop_thread_id = None
        self._slots = []

    def bind_loop_thread(self):
        """
        Must be called from inside the running event loop thread.
        Equivalent to Qt object thread affinity.
        """
        self._loop_thread_id = threading.get_ident()

    def connect(self, slot):
        if inspect.ismethod(slot):
            self._slots.append(weakref.WeakMethod(slot))
        else:
            self._slots.append(weakref.ref(slot))

    def emit(self, *args, **kwargs):
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

    def _invoke_direct(self, slot, *args, **kwargs):
        try:
            if inspect.iscoroutinefunction(slot):
                asyncio.create_task(slot(*args, **kwargs))
            else:
                slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (direct)")

    def _invoke_queued(self, slot, *args, **kwargs):
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
    def _safe_call(slot, *args, **kwargs):
        try:
            slot(*args, **kwargs)
        except Exception:
            logging.exception("Signal slot failed (queued)")
