import abc

from tavi.multithreading.proxy import Proxy


class TemplateModelInterface(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_next_file(self):
        pass


TemplateModelProxy = Proxy(TemplateModelInterface)
