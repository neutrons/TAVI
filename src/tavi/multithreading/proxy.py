import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Type, TypeVar

T = TypeVar("T")


def Proxy(_type: Type[T]):
    def __init__(self, host: T):
        self.host = host

    abstract_methods = getattr(_type, "__abstractmethods__", set())

    namespace = {"__init__": __init__}

    namespace = {"__init__": __init__}

    def make_proxy_method(method_name: str):
        def executeOnWorker(self, *args, **kwargs):
            with ThreadPoolExecutor(max_workers=2) as executor:
                host_method = getattr(self.host, method_name)
                fut = executor.submit(host_method, *args, **kwargs)
                return fut.result()

        # def executeOnWorker(self, *args, **kwargs):
        #     print(f"Running {method_name} on workerthread")
        #     return getattr(self.host, method_name)(*args, **kwargs)

        executeOnWorker.__name__ = method_name
        return executeOnWorker

    for name in abstract_methods:
        namespace[name] = make_proxy_method(name)

    ProxyClass = type(f"Proxy{_type.__name__}", (_type,), namespace)
    return ProxyClass