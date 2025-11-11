import abc

from tavi.multithreading.proxy import Proxy


class RandomModelInterface(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_next_file(self):
        pass


RandomModelProxy = Proxy(RandomModelInterface)
