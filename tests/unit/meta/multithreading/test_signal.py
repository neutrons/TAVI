import asyncio
import threading
import time
import unittest

from tavi.meta.multithreading.signal import Signal


def run_loop_once(loop):
    loop.call_soon(loop.stop)
    loop.run_forever()


def run_loop(loop, duration=0.05):
    loop.call_later(duration, loop.stop)
    loop.run_forever()


class TestSignal(unittest.TestCase):
    def setUp(self):
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()
        asyncio.set_event_loop(None)

    def test_direct_emit_executes_immediately(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = []

        def slot(value):
            called.append((value, threading.get_ident()))

        sig.connect(slot)
        sig.emit(123)

        self.assertEqual(called, [(123, threading.get_ident())])

    def test_queued_emit_from_worker_thread(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = []

        def slot(value):
            called.append((value, threading.get_ident()))

        sig.connect(slot)

        def worker():
            sig.emit(456)

        t = threading.Thread(target=worker)
        t.start()
        t.join()

        self.assertEqual(called, [])

        run_loop(self.loop)

        self.assertEqual(len(called), 1)
        value, thread_id = called[0]
        self.assertEqual(value, 456)
        self.assertEqual(thread_id, threading.get_ident())

    def test_async_slot_direct_requires_running_loop(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = asyncio.Event()

        async def slot(value):
            self.assertEqual(value, 1)
            called.set()

        sig.connect(slot)

        # Schedule emit while loop is running
        self.loop.call_soon(sig.emit, 1)
        run_loop(self.loop)

        self.loop.run_until_complete(
            asyncio.wait_for(called.wait(), timeout=0.5)
        )

    def test_async_slot_from_worker_thread(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = asyncio.Event()

        async def slot(value):
            called.set()

        sig.connect(slot)

        def worker():
            sig.emit(2)

        t = threading.Thread(target=worker)
        t.start()
        t.join()

        run_loop(self.loop)

        self.loop.run_until_complete(
            asyncio.wait_for(called.wait(), timeout=0.5)
        )

    def test_weakref_cleanup(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = []

        def make_slot():
            def slot():
                called.append(1)
            return slot

        slot = make_slot()
        sig.connect(slot)

        del slot

        sig.emit()
        run_loop(self.loop)

        self.assertEqual(called, [])
        self.assertEqual(sig._slots, [])

    def test_slot_order_preserved(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = []

        def slot1():
            called.append(1)

        def slot2():
            called.append(2)

        sig.connect(slot1)
        sig.connect(slot2)

        sig.emit()

        self.assertEqual(called, [1, 2])

    def test_slot_exception_does_not_break_signal(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        called = []

        def bad_slot():
            raise RuntimeError("boom")

        def good_slot():
            called.append("ok")

        sig.connect(bad_slot)
        sig.connect(good_slot)

        sig.emit()

        self.assertEqual(called, ["ok"])
        
    def test_queued_slot_executes_on_loop_thread(self):
        sig = Signal(self.loop)
        sig.bind_loop_thread()

        executed_thread_ids = []

        def slot():
            executed_thread_ids.append(threading.get_ident())

        sig.connect(slot)

        loop_thread_id = threading.get_ident()

        def worker():
            # Emit from a non-loop thread
            sig.emit()

        t = threading.Thread(target=worker)
        t.start()
        t.join()

        # Slot should not run until the loop processes events
        self.assertEqual(executed_thread_ids, [])

        run_loop(self.loop)

        self.assertEqual(len(executed_thread_ids), 1)
        self.assertEqual(executed_thread_ids[0], loop_thread_id)
