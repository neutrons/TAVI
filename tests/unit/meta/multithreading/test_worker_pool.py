import time
import threading
import unittest

from tavi.library.data.model_response import ModelResponse, ResponseCode


# -------------------------
# Minimal Signal test double
# -------------------------


class DummySignal:
    """
    Minimal replacement for Signal that synchronously
    executes connected callbacks.
    """

    def __init__(self, loop=None):
        self._callbacks = []
        self.emitted = False
        self.loopBound = False
        self.args = None

    def bind_loop_thread(self):
        self.loopBound = True

    def connect(self, callback):
        self._callbacks.append(callback)

    def emit(self, *args):
        self.emitted = True
        self.args = args
        for cb in self._callbacks:
            cb()


# -------------------------
# TestCase
# -------------------------


class TestWorkerAndWorkerPool(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Patch Signal once for the entire test class,
        then import the module under test.
        """
        import importlib
        import tavi.meta.multithreading.signal

        cls._original_signal = tavi.meta.multithreading.signal.Signal
        tavi.meta.multithreading.signal.Signal = DummySignal

        cls.worker_module = importlib.import_module("tavi.meta.multithreading.worker_pool")  # adjust if needed
        importlib.reload(cls.worker_module)

    @classmethod
    def tearDownClass(cls):
        """
        Restore original Signal implementation.
        """
        import tavi.meta.multithreading.signal

        tavi.meta.multithreading.signal.Signal = cls._original_signal

    def setUp(self):
        """
        Ensure WorkerPool is reset between tests.
        """
        WorkerPool = self.worker_module.WorkerPool
        pool = WorkerPool()
        pool.threads.clear()
        pool.worker_queue.clear()

    # -------------------------
    # Worker tests
    # -------------------------

    def test_worker_success_path(self):
        WorkerPool = self.worker_module.WorkerPool
        ResponseCode = self.worker_module.ResponseCode
        ModelResponse = self.worker_module.ModelResponse

        def target(a, b):
            return ModelResponse(code=ResponseCode.OK, data = a + b)

        worker = WorkerPool().create_worker(target, 2, 3)
        worker.run()

        self.assertTrue(worker.finished.loopBound)
        self.assertTrue(worker.finished.emitted)

    def test_worker_exception_path(self):
        WorkerPool = self.worker_module.WorkerPool

        def target():
            raise RuntimeError("failure")

        worker = WorkerPool().create_worker(target)
        worker.run()

        self.assertTrue(worker.finished.emitted)

    # -------------------------
    # WorkerPool tests
    # -------------------------

    def test_workerpool_singleton(self):
        WorkerPool = self.worker_module.WorkerPool

        pool1 = WorkerPool()
        pool2 = WorkerPool()

        self.assertIs(pool1, pool2)

    def test_workerpool_starts_thread(self):
        WorkerPool = self.worker_module.WorkerPool

        pool = WorkerPool()
        pool.max_threads = 1

        ran = threading.Event()

        def target():
            ran.set()
            time.sleep(0.05)
            return ModelResponse(code=ResponseCode.OK)

        worker = pool.create_worker(target)
        pool.submit_worker(worker)

        self.assertEqual(len(pool.threads), 1)

        ran.wait(timeout=1)
        time.sleep(0.05)

        self.assertEqual(len(pool.threads), 0)

    def test_workerpool_queues_and_dequeues(self):
        WorkerPool = self.worker_module.WorkerPool

        pool = WorkerPool()
        pool.max_threads = 1

        block = threading.Event()

        def slow_target():
            block.wait()
            return ModelResponse(code=ResponseCode.OK)

        def fast_target():
            return ModelResponse(code=ResponseCode.OK)

        w1 = pool.create_worker(slow_target)
        w2 = pool.create_worker(fast_target)

        pool.submit_worker(w1)
        pool.submit_worker(w2)

        self.assertEqual(len(pool.threads), 1)
        self.assertEqual(len(pool.worker_queue), 1)

        block.set()
        time.sleep(0.1)

        self.assertEqual(len(pool.worker_queue), 0)
        self.assertEqual(len(pool.threads), 0)

    def test_workerpool_finished_signal(self):
        WorkerPool = self.worker_module.WorkerPool
        ResponseCode = self.worker_module.ResponseCode
        ModelResponse = self.worker_module.ModelResponse

        def target(a, b):
            return ModelResponse(code=ResponseCode.OK, data = a + b)
        
        finishedCalled = False
        def slot():
            nonlocal finishedCalled
            finishedCalled = True

        worker = WorkerPool().create_worker(target, 2, 3)
        worker.finished.connect(slot)
        worker.run()

        self.assertTrue(worker.finished.emitted)

        assert finishedCalled