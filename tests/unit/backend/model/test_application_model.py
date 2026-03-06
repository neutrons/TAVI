

import pytest
import unittest
from unittest import mock, TestCase

from tavi.backend.model.application_model import ApplicationModel


class TestApplicationModel(TestCase):
    def setUp(self):
        self.mock_filestore = mock.Mock()
        self.application_model = ApplicationModel(self.mock_filestore)
        pass

    def tearDown(self):
        pass
    
    
    def test_write_error_log(self):
        self.application_model.write_error_log("test")
        
        self.mock_filestore.write_user_data_file.assert_called_once_with(mock.ANY, "test")