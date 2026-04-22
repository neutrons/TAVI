"""Unit tests for individual classification rules."""

import unittest
import unittest.mock

from neutrons_standard.config import Resource

from tavi.backend.classification.rule.ORNL.def_xy_rule import DefXYRule
from tavi.backend.classification.rule.ORNL.inst_in_filename_rule import (
    InstrumentInFilenameRule,
)
from tavi.backend.classification.rule.ORNL.spice_file_ext_rule import SpiceFileExtRule
from tavi.backend.classification.rule.dat_format_rule import DatFormatRule
from tavi.backend.classification.rule.hashtag_comment_rule import HashtagCommentRule
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class FakeFileStore(FileStoreInterface):
    """Fake implementation of FileStoreInterface for testing."""

    def __init__(self):
        pass

    def fetch_files_at(self, path: str) -> list[str]:
        pass

    def validate_file(self, file_path: str) -> bool:
        pass

    def write_user_data_file(self, subpath: str, value: str) -> None:
        """Not implemented for testing."""
        pass

    def write_text_file(self, path: str, value: str) -> None:
        """Not implemented for testing."""
        pass

    def read_text_file(self, path: str) -> str:
        """Read text from file."""
        with open(path, "r") as f:
            return f.read()

    def get_file_ext(self, path: str) -> str:
        """Get file extension from path."""
        return f".{path.split(".")[-1]}"

    def get_file_name(self, path: str) -> str:
        """Get file name from path."""
        return path.split("/")[-1]
    
    def get_file_size_mb(self, file_path: str) -> float:
        pass

    def get_parent(self, file_path: str) -> str:
        return "dummy"

    def join_path(self, root_path: str, target_path: str) -> str:
        return "dummy"


class TestInstrumentInFilenameRule(unittest.TestCase):
    """Test cases for InstrumentInFilenameRule."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_file_path = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.filestore = FakeFileStore()
        self.rule = InstrumentInFilenameRule()

    def test_hb1a_file_has_instrument_name(self) -> None:
        """Test that HB1A filename contains valid instrument name."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertEqual(
            score,
            1,
            "HB1A file should match instrument name rule",
        )
    
    def test_prints_debug(self) -> None:
        with unittest.mock.patch("tavi.backend.classification.rule.ORNL.inst_in_filename_rule.logger") as mock_logger:
            score = self.rule.get_score(None, None)
            self.assertEqual(
                score,
                0,
                "Bad input should score a 0.",
            )
            mock_logger.debug.assert_called()
        
    def test_rule_returns_integer(self) -> None:
        """Test that rule returns an integer score."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertIsInstance(score, int, "Rule should return an integer score")


class TestSpiceFileExtRule(unittest.TestCase):
    """Test cases for SpiceFileExtRule."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_file_path = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.filestore = FakeFileStore()
        self.rule = SpiceFileExtRule()

    def test_dat_file_extension_match(self) -> None:
        """Test that .dat file extension is recognized."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertEqual(
            score,
            1,
            "File with .dat extension should match SPICE file extension rule",
        )
        
    def test_prints_debug(self) -> None:
        with unittest.mock.patch("tavi.backend.classification.rule.ORNL.spice_file_ext_rule.logger") as mock_logger:
            score = self.rule.get_score(None, None)
            self.assertEqual(
                score,
                0,
                "Bad input should score a 0.",
            )
            mock_logger.debug.assert_called()

    def test_rule_returns_integer(self) -> None:
        """Test that rule returns an integer score."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertIsInstance(score, int, "Rule should return an integer score")


class TestDatFormatRule(unittest.TestCase):
    """Test cases for DatFormatRule."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_file_path_HB1A = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.test_file_path_HBZ = Resource.getPath(
            "inputs/rule_based_classifier/033947"
        )
        self.test_file_path_garbage = Resource.getPath(
            "inputs/rule_based_classifier/garbage"
        )
        self.filestore = FakeFileStore()
        self.rule = DatFormatRule()

    def test_hb1a_file_is_dat_format(self) -> None:
        """Test that HB1A file is in DAT format."""
        score = self.rule.get_score(self.test_file_path_HB1A, self.filestore)
        self.assertEqual(
            score,
            1,
            "HB1A file should match DAT format rule",
        )

    def test_fail(self) -> None:
        """Test that a garbage file fails this rule."""
        score = self.rule.get_score(self.test_file_path_garbage, self.filestore)
        self.assertEqual(
            score,
            0,
            "Garbage file should not match DAT format rule",
        )

    def test_hbz_file_is_dat_format(self) -> None:
        """Test that HB1A file is in DAT format."""
        score = self.rule.get_score(self.test_file_path_HBZ, self.filestore)
        self.assertEqual(
            score,
            1,
            "HBZ file should match DAT format rule",
        )

    def test_rule_returns_integer(self) -> None:
        """Test that rule returns an integer score."""
        score = self.rule.get_score(self.test_file_path_HB1A, self.filestore)
        self.assertIsInstance(score, int, "Rule should return an integer score")


class TestHashtagCommentRule(unittest.TestCase):
    """Test cases for HashtagCommentRule."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_file_path = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.filestore = FakeFileStore()
        self.rule = HashtagCommentRule()

    def test_hb1a_file_has_hashtag_comments(self) -> None:
        """Test that HB1A file contains hashtag comments."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertEqual(
            score,
            1,
            "HB1A file should match hashtag comment rule",
        )

    def test_prints_debug(self) -> None:
        with unittest.mock.patch("tavi.backend.classification.rule.hashtag_comment_rule.logger") as mock_logger:
            score = self.rule.get_score(None, None)
            self.assertEqual(
                score,
                0,
                "Bad input should score a 0.",
            )
            mock_logger.debug.assert_called()

    def test_rule_returns_integer(self) -> None:
        """Test that rule returns an integer score."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertIsInstance(score, int, "Rule should return an integer score")


class TestDefXYRule(unittest.TestCase):
    """Test cases for DefXYRule."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.test_file_path = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.filestore = FakeFileStore()
        self.rule = DefXYRule()

    def test_hb1a_file_has_def_metadata_fields(self) -> None:
        """Test that HB1A file contains DEF metadata fields."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertEqual(
            score,
            1,
            "HB1A file should match DEF metadata field rule",
        )

    def test_rule_returns_integer(self) -> None:
        """Test that rule returns an integer score."""
        score = self.rule.get_score(self.test_file_path, self.filestore)
        self.assertIsInstance(score, int, "Rule should return an integer score")


if __name__ == "__main__":
    unittest.main()
