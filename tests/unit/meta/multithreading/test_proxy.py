import abc
from unittest import mock, TestCase

import pytest

from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.meta.multithreading.proxy import Proxy


class SampleAbstractClass(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def some_method(self):
        pass


class SampleAbstractClassImpl(SampleAbstractClass):
    test_var = lambda: True

    def some_method(self):
        self.test_var()
        return ModelResponse(code=ResponseCode.OK)

    def some_non_abstract_method(self):
        pass


SampleAbstractClassProxy = Proxy(SampleAbstractClass)


class SampleSyncAbstractClass(metaclass=abc.ABCMeta):
    _proxy_sync_methods = frozenset({"sync_method"})

    @abc.abstractmethod
    def sync_method(self, value):
        pass

    @abc.abstractmethod
    def async_method(self):
        pass


class SampleSyncAbstractClassImpl(SampleSyncAbstractClass):
    def sync_method(self, value):
        return value * 2

    def async_method(self):
        return ModelResponse(code=ResponseCode.OK)


SampleSyncAbstractClassProxy = Proxy(SampleSyncAbstractClass)


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
            assert self.inst.worker_pool is mock_worker_pool
            assert mock_worker_pool.create_worker.called
            assert mock_worker_pool.submit_worker.called

    def test_non_interface_method(self):
        with pytest.raises(AttributeError, match="object has no attribute"):
            self.inst.some_non_abstract_method()
        self.inst.host.some_non_abstract_method()


class TestSyncProxyMethod(TestCase):
    def setUp(self):
        self.inst = SampleSyncAbstractClassProxy(SampleSyncAbstractClassImpl())

    def test_sync_method_returns_value_directly(self):
        assert self.inst.sync_method(21) == 42

    def test_sync_method_does_not_use_worker_pool(self):
        with mock.patch.object(self.inst, "worker_pool", mock.Mock()) as mock_worker_pool:
            self.inst.sync_method(1)
            assert not mock_worker_pool.create_worker.called
            assert not mock_worker_pool.submit_worker.called

    def test_async_method_still_uses_worker_pool(self):
        with mock.patch.object(self.inst, "worker_pool", mock.Mock()) as mock_worker_pool:
            self.inst.async_method()
            assert mock_worker_pool.create_worker.called
            assert mock_worker_pool.submit_worker.called
