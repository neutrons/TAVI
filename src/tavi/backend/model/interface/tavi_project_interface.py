import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.meta.multithreading.proxy import Proxy


class TaviProjectInterface(Model, metaclass=abc.ABCMeta):
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
