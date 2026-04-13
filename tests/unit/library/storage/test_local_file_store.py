

from pathlib import Path
from tempfile import TemporaryDirectory
import pytest
import unittest
from unittest import mock, TestCase


from tavi.library.storage.local_file_store import LocalFileStore
from util.Config_helpers import Config_override
from neutrons_standard.config import Config


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

    def test_read_text_file(self):
        with TemporaryDirectory() as tmp:
            filename = "test_read.txt"
            filepath = Path(f"{tmp}/{filename}")
            test_content = "Hello, World!"
            self.local_filestore.write_text_file(str(filepath), test_content)
            
            read_content = self.local_filestore.read_text_file(str(filepath))
            assert read_content == test_content

    def test_read_text_file_raises_for_nonexistent_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/nonexistent.txt")
            
            with pytest.raises(RuntimeError):
                self.local_filestore.read_text_file(str(filepath))

    def test_get_file_ext(self):
        with TemporaryDirectory() as tmp:
            filename = "test_file.txt"
            filepath = Path(f"{tmp}/{filename}")
            self.local_filestore.write_text_file(str(filepath), "content")
            
            ext = self.local_filestore.get_file_ext(str(filepath))
            assert ext == ".txt"

    def test_get_file_ext_various_extensions(self):
        with TemporaryDirectory() as tmp:
            test_cases = [
                ("file.py", ".py"),
                ("data.json", ".json"),
                ("document.pdf", ".pdf"),
                ("archive.tar.gz", ".gz"),
            ]
            
            for filename, expected_ext in test_cases:
                filepath = Path(f"{tmp}/{filename}")
                self.local_filestore.write_text_file(str(filepath), "content")
                
                ext = self.local_filestore.get_file_ext(str(filepath))
                assert ext == expected_ext

    def test_get_file_ext_raises_for_nonexistent_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/nonexistent.txt")
            
            with pytest.raises(RuntimeError):
                self.local_filestore.get_file_ext(str(filepath))

    def test_get_file_size_mb(self):
        with TemporaryDirectory() as tmp:
            filename = "test_file.txt"
            filepath = Path(f"{tmp}/{filename}")
            test_content = "x" * 1024  # 1 KB
            self.local_filestore.write_text_file(str(filepath), test_content)
            
            size_mb = self.local_filestore.get_file_size_mb(str(filepath))
            expected_mb = 1024 / (1024 * 1024)
            assert abs(size_mb - expected_mb) < 0.001

    def test_get_file_size_mb_raises_for_nonexistent_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/nonexistent.txt")
            
            with pytest.raises(RuntimeError):
                self.local_filestore.get_file_size_mb(str(filepath))

    def test_validate_file_with_valid_file(self):
        with TemporaryDirectory() as tmp:
            with Config_override("library.filestore.raw.size-limit", 100):
                filename = "valid_file.txt"
                filepath = Path(f"{tmp}/{filename}")
                self.local_filestore.write_text_file(str(filepath), "content")
                
                is_valid = self.local_filestore.validate_file(str(filepath))
                assert is_valid is True

    def test_validate_file_with_nonexistent_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/nonexistent.txt")
            
            is_valid = self.local_filestore.validate_file(str(filepath))
            assert is_valid is False

    def test_validate_file_with_directory(self):
        with TemporaryDirectory() as tmp:
            is_valid = self.local_filestore.validate_file(tmp)
            assert is_valid is False

    def test_validate_file_respects_size_limit(self):
        with TemporaryDirectory() as tmp:
            filename = "large_file.txt"
            filepath = Path(f"{tmp}/{filename}")
            
            # Create a file that exceeds the size limit
            size_limit_mb = 0.001  # 1 KB limit
            with Config_override("library.filestore.raw.size-limit", size_limit_mb):
                size_bytes = int(size_limit_mb * 1024 * 1024) + 1024  # 1 KB over limit
                self.local_filestore.write_text_file(str(filepath), "x" * size_bytes)
                
                is_valid = self.local_filestore.validate_file(str(filepath))
                assert is_valid is False

    def test_fetch_files_at_with_valid_directory(self):
        with TemporaryDirectory() as tmp:
            with Config_override("library.filestore.raw.size-limit", 100):
                # Create multiple test files
                files = ["file1.txt", "file2.txt", "file3.txt"]
                for filename in files:
                    filepath = Path(f"{tmp}/{filename}")
                    self.local_filestore.write_text_file(str(filepath), f"content of {filename}")
                
                fetched_files = self.local_filestore.fetch_files_at(tmp)
                assert len(fetched_files) == 3
                for filepath in fetched_files:
                    assert Path(filepath).name in files

    def test_fetch_files_at_excludes_nonexistent_files(self):
        with TemporaryDirectory() as tmp:
            with Config_override("library.filestore.raw.size-limit", 100):
                # Create a valid file
                filepath1 = Path(f"{tmp}/valid.txt")
                self.local_filestore.write_text_file(str(filepath1), "content")
                
                # Create a directory (invalid)
                subdir = Path(f"{tmp}/subdir")
                subdir.mkdir()
                
                fetched_files = self.local_filestore.fetch_files_at(tmp)
                # Should only include the valid file
                assert len(fetched_files) == 1
                assert Path(fetched_files[0]).name == "valid.txt"

    def test_fetch_files_at_raises_for_nonexistent_path(self):
        nonexistent_path = "/tmp/nonexistent_path_12345"
        
        with pytest.raises(RuntimeError, match="does not exist"):
            self.local_filestore.fetch_files_at(nonexistent_path)

    def test_fetch_files_at_raises_for_file_path(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/test.txt")
            self.local_filestore.write_text_file(str(filepath), "content")
            
            with pytest.raises(RuntimeError, match="is a file"):
                self.local_filestore.fetch_files_at(str(filepath))

    def test_write_text_file_overwrites_existing_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/test.txt")
            self.local_filestore.write_text_file(str(filepath), "first content")
            self.local_filestore.write_text_file(str(filepath), "second content")
            
            read_content = self.local_filestore.read_text_file(str(filepath))
            assert read_content == "second content"

    def test_write_user_data_file_rejects_value_with_leading_slash(self):
        with TemporaryDirectory() as tmp:
            with Config_override("user.application.home", tmp):
                # Note: The actual code checks if value (content) starts with "/" not the path
                with pytest.raises(RuntimeError, match="should not start with a /"):
                    self.local_filestore.write_user_data_file("file.txt", "/content")

    def test_is_real_file_returns_true_for_valid_file(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(f"{tmp}/test.txt")
            self.local_filestore.write_text_file(str(filepath), "content")
            
            result = self.local_filestore._is_real_file(filepath)
            assert result is True

    def test_is_real_file_returns_false_for_nonexistent_file(self):
        filepath = Path("/tmp/nonexistent_file_12345.txt")
        
        result = self.local_filestore._is_real_file(filepath)
        assert result is False

    def test_is_real_file_returns_false_for_directory(self):
        with TemporaryDirectory() as tmp:
            filepath = Path(tmp)
            
            result = self.local_filestore._is_real_file(filepath)
            assert result is False

    def test_is_real_file_throws_when_requested_for_invalid_file(self):
        filepath = Path("/tmp/nonexistent_file_12345.txt")
        
        with pytest.raises(RuntimeError, match="does not exist or is not a file"):
            self.local_filestore._is_real_file(filepath, throws=True)