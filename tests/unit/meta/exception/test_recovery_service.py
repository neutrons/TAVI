

from pathlib import Path
from tempfile import TemporaryDirectory
import pytest
import unittest
from unittest import mock, TestCase


from tavi.library.storage.local_file_store import LocalFileStore
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.exception.recovery_service import RecoveryService
from tavi.meta.exception.tavi_exception import TaviError
from util.Config_helpers import Config_override


class DummyError(TaviError):
    pass


class TestRecoveryService(TestCase):
    def setUp(self):
        self.recovery_service = RecoveryService()

    def tearDown(self):
        pass
    
    def test_register(self):
        
        def handle(_:DummyError):
            pass
        
        self.recovery_service.register(TaviError, handle)
        
        assert len(self.recovery_service.exception_handlers) == 1
        
        with pytest.raises(RuntimeError, match=f"Only Exceptions of subtype {TaviError.__name__} can be registered."):
            self.recovery_service.register(ValueError, handle)

    def test_handle_exception(self):
        didHandle = False
        
        def handle(ex: DummyError):
            nonlocal didHandle
            didHandle = True
            assert isinstance(ex, DummyError)
            
        self.recovery_service.register(DummyError, handle)
        
        self.recovery_service.handle_exception(ExceptionEvent(error=DummyError("message", "stack_trace")))
        
        assert didHandle
        
    def test_default_handler(self):
        with pytest.raises(RuntimeError, match="FATAL: no handler for exception found"):
            self.recovery_service.handle_exception(ExceptionEvent(error=DummyError("message", "stack_trace")))