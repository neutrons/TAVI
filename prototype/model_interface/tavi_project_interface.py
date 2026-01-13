import abc

from tavi.multithreading.proxy import Proxy


class TaviProjectInterface(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def set_selected_scan(self):
        pass

    @abc.abstractmethod
    def get_selected_metadata(self):
        pass

    @abc.abstractmethod
    def load_manager(self):
        pass

    @abc.abstractmethod
    def load(self):
        pass


TaviProjectProxy = Proxy(TaviProjectInterface)
