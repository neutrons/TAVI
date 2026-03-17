import abc
from typing import Any
from unittest import mock, TestCase

import pytest


from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.data.scan import ScanMetadata
from tavi.library.data.scan import Scan
from tavi.library.data.scan import ScanData
from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.loader_registry import LoaderRegistry



class DummyLoader(AbstractLoader):
    
    def __init__(self, filestore:Filestore):
        super().__init__(filestore)
    
    def load(path:str) -> Scan:
        pass
    
    def get_scan_type(self) -> RawScanType:
        """Get the scan type identifier for this loader."""
        pass
    
    def get_score(path:str) -> int:
        pass
    
    def parse_metadata(path:str) -> ScanMetadata:
        pass
    
    def parse_scan_values(path:str) -> ScanData:
        pass
    
    def parse_external_metadata(path:str) -> dict[str, Any]:
        pass
    
    def adapt_scan_data(meta: ScanMetadata, values:ScanData):
        pass

class SmartyLoader(DummyLoader):
    pass


class TestFileOperations(TestCase):
    def setUp(self):
        self.inst: LoaderRegistry = LoaderRegistry(None)

    def tearDown(self):
        pass

    def test_register(self):
        mock_filestore_1 = mock.Mock()
        
        loader = DummyLoader(mock_filestore_1)
        
        self.inst.register(DummyLoader.__name__, loader)
        
        assert self.inst.registry[DummyLoader.__name__] is loader
        
    def test_get_loaders(self):
        mock_filestore_1 = mock.Mock()
        mock_filestore_2 = mock.Mock()
        
        loader_1 = DummyLoader(mock_filestore_1)
        loader_2 = SmartyLoader(mock_filestore_2)
        
        self.inst.set_filestore(mock_filestore_1)
        
        self.inst.register(DummyLoader.__name__, loader_1)
        self.inst.register(SmartyLoader.__name__, loader_2)
        
        assert loader_1 in self.inst.get_loaders()
        assert loader_2 in self.inst.get_loaders()
        
        for loader in self.inst.get_loaders():
            assert loader.filestore is mock_filestore_1
            
        self.inst.set_filestore(mock_filestore_2)
        
        for loader in self.inst.get_loaders():
            assert loader.filestore is mock_filestore_2