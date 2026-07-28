"""Proxy."""

from typing import Any, Callable, Type, TypeVar

from tavi.meta.multithreading.worker_pool import WorkerPool

# include this in a config file
num_of_workers = 1

T = TypeVar("T")


def Proxy(_type: Type[T]) -> type:
    """Generate a wrapper class that runs method calls to host on a thread."""

    def __init__(self: Any, host: T) -> None:
        self.host = host
        self.worker_pool = WorkerPool()

    abstract_methods = getattr(_type, "__abstractmethods__", set())
    # Methods that are pure, side-effect-free reads of already-loaded data may opt out of the
    # async worker-thread dispatch above by listing their name here, on the interface class.
    sync_methods = getattr(_type, "_proxy_sync_methods", frozenset())

    namespace = {"__init__": __init__}

    def make_proxy_method(method_name: str) -> Callable:
        def executeOnWorker(self: T, *args: Any, **kwargs: Any) -> None:
            host_method = getattr(self.host, method_name)
            worker = self.worker_pool.create_worker(host_method, *args, **kwargs)
            self.worker_pool.submit_worker(worker)
            # TODO: move to debug logger
            print("sent to worker!")
            # Models do not return values.
            # Returning values is synchronous behavior that will freeze the main thread.
            # They may submit updates via events.
            # This also helps ensure *all* consumers are synchronized.

        executeOnWorker.__name__ = method_name
        return executeOnWorker

    def make_sync_proxy_method(method_name: str) -> Callable:
        def callDirectly(self: T, *args: Any, **kwargs: Any) -> Any:
            return getattr(self.host, method_name)(*args, **kwargs)

        callDirectly.__name__ = method_name
        return callDirectly

    for name in abstract_methods:
        namespace[name] = make_sync_proxy_method(name) if name in sync_methods else make_proxy_method(name)

    ProxyClass = type(f"Proxy{_type.__name__}", (_type,), namespace)
    return ProxyClass
