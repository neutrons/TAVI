import abc
from unittest import mock, TestCase

import pytest

from tavi.meta.multithreading.proxy import Proxy


class SampleAbstractClass(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def some_method(self):
        pass


class SampleAbstractClassImpl(SampleAbstractClass):
    test_var = lambda: True

    def some_method(self):
        self.test_var()

    def some_non_abstract_method(self):
        pass


SampleAbstractClassProxy = Proxy(SampleAbstractClass)


class TestFileOperations(TestCase):
    def setUp(self):
        self.inst = SampleAbstractClassProxy(SampleAbstractClassImpl())

    def tearDown(self):
        pass

    def test_init(self):
        assert isinstance(self.inst, SampleAbstractClass)
        assert not isinstance(self.inst, SampleAbstractClassImpl)

        assert isinstance(self.inst.host, SampleAbstractClass)
        assert isinstance(self.inst.host, SampleAbstractClassImpl)

    def test_method_exec(self):
        self.inst.host.test_var = mock.Mock()
        with mock.patch.object(self.inst.host, "test_var", mock.Mock()) as test_var_mock:
            self.inst.some_method()
            assert test_var_mock.called

    def test_method_worker_delegation(self):
        with mock.patch.object(self.inst, "worker_pool", mock.Mock()) as mock_worker_pool:
            self.inst.some_method()
            assert mock_worker_pool.create_worker.called
            assert mock_worker_pool.submit_worker.called

    def test_non_interface_method(self):
        with pytest.raises(AttributeError, match="object has no attribute"):
            self.inst.some_non_abstract_method()
        self.inst.host.some_non_abstract_method()
