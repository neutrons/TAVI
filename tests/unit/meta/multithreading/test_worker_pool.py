import time
import threading
import unittest


# -------------------------
# Minimal Signal test double
# -------------------------

class TestSignal:
    """
    Minimal replacement for Signal that synchronously
    executes connected callbacks.
    """
    def __init__(self, loop=None):
        self._callbacks = []
        self.emitted = False
        self.args = None

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
        tavi.meta.multithreading.signal.Signal = TestSignal

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
        Worker = self.worker_module.Worker
        ResponseCode = self.worker_module.ResponseCode
        ModelResponse = self.worker_module.ModelResponse

        def target(a, b):
            return a + b

        worker = Worker(target, 2, 3)
        worker.run()

        self.assertTrue(worker.result.emitted)
        self.assertTrue(worker.success.emitted)
        self.assertTrue(worker.finished.emitted)

        result = worker.result.args[0]
        self.assertIsInstance(result, ModelResponse)
        self.assertEqual(result.code, ResponseCode.OK)
        self.assertEqual(result.data, 5)

    def test_worker_exception_path(self):
        Worker = self.worker_module.Worker
        ResponseCode = self.worker_module.ResponseCode
        ModelResponse = self.worker_module.ModelResponse

        def target():
            raise RuntimeError("failure")

        worker = Worker(target)
        worker.run()

        self.assertTrue(worker.result.emitted)
        self.assertTrue(worker.success.emitted)
        self.assertTrue(worker.finished.emitted)

        result = worker.result.args[0]
        self.assertIsInstance(result, ModelResponse)
        self.assertEqual(result.code, ResponseCode.ERROR)
        self.assertIn("failure", result.message)

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
            return "ok"

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
            return "slow"

        def fast_target():
            return "fast"

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

