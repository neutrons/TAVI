

from pathlib import Path
from tempfile import TemporaryDirectory
import pytest
import unittest
from unittest import mock, TestCase


from tavi.library.storage.local_file_store import LocalFileStore
from tests.util.Config_helpers import Config_override


class TestLocalFilestore(TestCase):
    def setUp(self):
        self.local_filestore = LocalFileStore()

    def tearDown(self):
        pass
    
    
    def test_write_text_file(self):
        with TemporaryDirectory() as tmp:
                filename = "delete_me.test"
                filepath = Path(f"{tmp}/{filename}")
                self.local_filestore.write_text_file(str(filepath), "test file.")
                assert filepath.exists()
    
    def test_write_user_data_file(self):
        with TemporaryDirectory() as tmp:
            with Config_override("user.application.home", tmp):
                filename = "delete_me.test"
                self.local_filestore.write_user_data_file(filename, "test file.")
                filepath = Path(f"{tmp}/{filename}")
                assert filepath.exists()