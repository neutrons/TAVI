from concurrent.futures import ThreadPoolExecutor
from typing import Type, TypeVar

from tavi.meta.multithreading.worker_pool import WorkerPool

# include this in a config file
num_of_workers = 1

T = TypeVar("T")

def Proxy(_type: Type[T]):
    def __init__(self, host: T):
        self.host = host
        self.worker_pool = WorkerPool()

    abstract_methods = getattr(_type, "__abstractmethods__", set())

    namespace = {"__init__": __init__}

    def make_proxy_method(method_name: str):
        def executeOnWorker(self, *args, **kwargs):
            host_method = getattr(self.host, method_name)
            worker = self.worker_pool.create_worker(host_method, *args, **kwargs)
            self.worker_pool.submit_worker(worker)
            # TODO: move to debug logger
            print("sent to worker!")
            # Models do not return values.
            # Returning values is synchonous behavior that will freeze the main thread.
            # They may submit updates via events.
            # This also helps ensure *all* consumers are synchronized.

        executeOnWorker.__name__ = method_name
        return executeOnWorker

    for name in abstract_methods:
        namespace[name] = make_proxy_method(name)

    ProxyClass = type(f"Proxy{_type.__name__}", (_type,), namespace)
    return ProxyClass
