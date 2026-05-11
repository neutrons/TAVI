"""Unit tests for RuleBasedClassifier and ORNLSpiceRuleSet."""

import unittest
from unittest import mock

from neutrons_standard.config import Resource

from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class TestFileStore(FileStoreInterface):
    """Test implementation of FileStoreInterface that reads from actual files."""

    def fetch_files_at(self, path: str) -> list[str]:
        pass

    def validate_file(self, file_path: str) -> bool:
        pass

    def write_user_data_file(self, subpath: str, value: str) -> None:
        """Not implemented for testing."""
        pass
    
    def read_user_data_file(self, file_subpath):
        return super().read_user_data_file(file_subpath)

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
        pass

    def join_path(self, root_path: str, target_path: str) -> str:
        pass


class TestRuleBasedClassifier(unittest.TestCase):
    """Test cases for RuleBasedClassifier using ORNLSpiceRuleSet."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        # Get the test data file using neutrons_standard.config.Resource
        self.test_file_path = Resource.getPath(
            "inputs/rule_based_classifier/HB1A_exp0978_scan0001.dat"
        )
        self.filestore = TestFileStore()
        self.classifier = RuleBasedClassifier(self.filestore)
        self.rule_set = ORNLSpiceRuleSet()

    def test_classifier_with_hb1a_file(self) -> None:
        """Test classifier with HB1A SPICE file.

        The HB1A file should score well because it has:
        - Instrument name (HB1A) in filename
        - .dat file extension
        - Hashtag comments (# prefix)
        - DAT format with space-separated values
        - DEF metadata fields

        """
        score = self.classifier.get_score(str(self.test_file_path), self.rule_set)

        # The score should be positive since the file matches most rules
        self.assertGreater(
            score,
            0,
            "HB1A file should score positive with ORNL SPICE rule set",
        )

    def test_classifier_score_components(self) -> None:
        """Test that individual rules contribute to the overall score."""
        score = self.classifier.get_score(str(self.test_file_path), self.rule_set)

        # Each rule that passes should contribute its weight (1 point)
        # We expect at least 5 rules to pass: instrument name, dat format, hashtag
        self.assertGreaterEqual(
            score,
            1.0,
            "HB1A file should match at least 5 rules",
        )

    def test_file_extension_matching(self) -> None:
        """Test that file extension is correctly identified."""
        ext = self.filestore.get_file_ext(str(self.test_file_path))
        self.assertEqual(ext, ".dat", "HB1A file should have .dat extension")

    def test_file_name_extraction(self) -> None:
        """Test that file name is correctly extracted."""
        name = self.filestore.get_file_name(str(self.test_file_path))
        self.assertEqual(
            name,
            "HB1A_exp0978_scan0001.dat",
            "File name should be correctly extracted",
        )

    def test_file_contains_hashtag_comments(self) -> None:
        """Test that the HB1A file contains hashtag comments."""
        content = self.filestore.read_text_file(str(self.test_file_path))
        lines = content.split("\n")
        hashtag_lines = [line for line in lines if line.startswith("#")]
        self.assertGreater(
            len(hashtag_lines),
            0,
            "HB1A file should contain lines starting with #",
        )

    def test_rule_set_has_all_rules(self) -> None:
        """Test that ORNLSpiceRuleSet contains all expected rules."""
        rules = self.rule_set.get_rules()
        self.assertEqual(
            len(rules),
            5,
            "ORNLSpiceRuleSet should have 5 rules",
        )

    def test_rule_weights(self) -> None:
        """Test that all rules have the correct weight."""
        for rule in self.rule_set.get_rules():
            weight = self.rule_set.get_weight(rule)
            self.assertEqual(
                weight,
                .2,
                f"Rule {rule.__class__.__name__} should have weight 1",
            )


if __name__ == "__main__":
    unittest.main()
