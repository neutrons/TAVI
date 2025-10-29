from typing import TypeVar, Generic, Type
import abc


T = TypeVar('T')

def Proxy(_type: Type[T]):
    def __init__(self, host:T):
        self.host = host
    
    abstract_methods = getattr(_type, "__abstractmethods__", set())

    namespace = {"__init__" : __init__}
    def make_proxy_method(method_name: str):
        def executeOnWorker(self, *args, **kwargs):
            print(f"Running {method_name} on workerthread")
            return getattr(self.host, method_name)(*args, **kwargs)
        executeOnWorker.__name__ = method_name
        return executeOnWorker
    
    for name in abstract_methods:
        namespace[name] = make_proxy_method(name)
    
    
    ProxyClass = type(f"Proxy{_type.__name__}", (_type,), namespace)
    return ProxyClass

# same file
# class MyClassInterface(metaclass=abc.ABCMeta):
#         @abc.abstractmethod
#         def foo(self):
#             pass
# MyClassProxy = Proxy(MyClassInterface)

# separate file                
# class MyClass(MyClassInterface):
#     def foo(self):
#         print("Foo!")




# if __name__ == "__main__":
#     myClass = MyClass()
#     myClassProxy = MyClassProxy(myClass)

#     print(type(myClass))
#     myClass.foo()
#     print(type(myClassProxy))
#     myClassProxy.foo()
