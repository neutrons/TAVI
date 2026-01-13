import abc

from tavi.backend.model.interface.model_interface import Model
from tavi.meta.multithreading.proxy import Proxy


class TemplateModelInterface(Model, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_next_file(self):
        pass


TemplateModelProxy = Proxy(TemplateModelInterface)
