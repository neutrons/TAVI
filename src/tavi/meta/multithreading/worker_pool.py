"""Worker and Worker Pool."""

import asyncio
import traceback
from threading import Thread
from typing import Any, Callable, Dict, List

from neutrons_standard.decorators.singleton import Singleton

from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError
from tavi.meta.multithreading.signal import Signal

# logger = taviLogger.getLogger(__name__)


class Worker:
    """Wrapper for threaded tasks that emits signals and handles errors."""

    target = None
    args = None

    def __init__(
        self, loop: asyncio.AbstractEventLoop, target: Callable, call_stack: str, *args: Any, **kwargs: Any
    ) -> None:
        """Initialize worker and setup expected async signals."""
        super().__init__()
        self.target = target
        self.call_stack = call_stack
        self.args = args
        self.kwargs = kwargs
        self.event_broker = EventBroker()
        self.finished: Signal = Signal(loop)  # None

    def bindSignals(self) -> None:
        """Bind the signals for use."""
        self.finished.bind_loop_thread()

    def run(self) -> None:
        """Long-running task."""
        results = None
        self.bindSignals()
        try:
            # Expects the return to be wrapped in a ModelResponse
            results: ModelResponse = self.target(*self.args, **self.kwargs)
        except Exception as e:  # noqa: BLE001
            import traceback

            stack_trace = f"{self.call_stack} \n {traceback.format_exc()}"
            error_message = str(e)
            results = ModelResponse(code=ResponseCode.ERROR, message=error_message)
            self.event_broker.publish(ExceptionEvent(e=NonRecoverableError(error_message, stack_trace)))

        self.finished.emit()
        if not isinstance(results, ModelResponse):
            raise ValueError("Worker target must return a ModelResponse object.")


@Singleton
class WorkerPool:
    """Creates and manages threads to run Workers on."""

    max_threads = 8
    threads: Dict[Worker, Thread] = {}
    worker_queue: List[Worker] = []

    def __init__(self) -> None:
        """Set up event loop so that Signals may work correctly."""
        self.loop = asyncio.new_event_loop()

    def create_worker(self, target: Callable, *args: Any, **kwargs: Any) -> Worker:
        """Create a worker."""
        call_stack = "".join(traceback.format_stack())
        return Worker(self.loop, target, call_stack, *args, **kwargs)

    def _dequeue_worker(self, worker: Worker) -> None:
        """Dequeues worker and starts the next in queue if it exists."""
        self.threads.pop(worker)
        if len(self.worker_queue) > 0:
            self.submit_worker(self.worker_queue.pop())

    def submit_worker(self, worker: Worker) -> None:
        """Queues or submits worker to thread."""
        if len(self.threads) >= self.max_threads:
            # add to queue
            self.worker_queue.append(worker)
        else:
            # spawn thread and delegate
            thread = Thread(target=worker.run)
            self.threads[worker] = thread

            worker.finished.connect(lambda: self._dequeue_worker(worker))

            # Start the thread
            thread.start()
