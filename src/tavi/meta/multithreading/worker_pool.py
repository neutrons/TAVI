
import asyncio
from threading import Thread
from typing import Dict, List

from tavi.meta.multithreading.signal import Signal
from tavi.meta.decorators.Singleton import Singleton
from tavi.library.data.ModelResponse import ModelResponse, ResponseCode

# logger = taviLogger.getLogger(__name__)


class Worker:
    target = None
    args = None

    def __init__(self, target, *args, **kwargs):
        super().__init__()
        self.target = target
        self.args = args
        self.kwargs = kwargs
        loop = asyncio.get_event_loop()
        self.finished = Signal(loop) # None
        self.success = Signal(loop) # bool
        self.result = Signal(loop) # object
        self.progress = Signal(loop) # int

    def run(self):
        """Long-running task."""
        results = None
        try:
            ret = self.target(*self.args, **self.kwargs)
            results = ModelResponse(code=ResponseCode.OK, data=ret)
            # results.code = 200 # set to 200 for testing
        except Exception as e:  # noqa: BLE001
            # logger.error(e)
            # if logger.isEnabledFor(logging.DEBUG):
            #     # print stacktrace
            #     import traceback

            #     traceback.print_exc()

            results = ModelResponse(code=ResponseCode.ERROR, message=str(e))

        if isinstance(results, ModelResponse):
            self.result.emit(results)
            self.success.emit(results.code > ResponseCode.OK)
            self.finished.emit()
        else:
            self.finished.emit()
            raise ValueError("Worker target must return a ModelResponse object.")


@Singleton
class WorkerPool:
    max_threads = 8
    threads: Dict[Worker, Thread] = {}
    worker_queue: List[Worker] = []
    
    def __init__(self):
        import asyncio

        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        

    def createWorker(self, target, *args, **kwargs):
        return Worker(target, *args, **kwargs)

    def _dequeueWorker(self, worker):
        self.threads.pop(worker)
        if len(self.worker_queue) > 0:
            self.submitWorker(self.worker_queue.pop())

    def submitWorker(self, worker):
        if len(self.threads) >= self.max_threads:
            # add to queue
            self.worker_queue.append(worker)
        else:
            # spawn thread and delegate
            thread = Thread(target=worker.run)
            self.threads[worker] = thread

            worker.finished.connect(lambda: self._dequeueWorker(worker))
            
            # Start the thread
            thread.start()