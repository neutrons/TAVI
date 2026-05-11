"""Unit tests for RawScanClassifier."""

import unittest
from unittest import mock

from tavi.backend.classification.raw_scan_classifier import RawScanClassifier
from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.loader_registry import LoaderRegistry
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class MockFileStore(FileStoreInterface):
    """Mock implementation of FileStoreInterface for testing."""

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
        """Not implemented for testing."""
        pass

    def get_file_ext(self, path: str) -> str:
        """Get file extension from path."""
        return path.split(".")[-1]

    def get_file_name(self, path: str) -> str:
        """Get file name from path."""
        return path.split("/")[-1]
    
    def get_file_size_mb(self, file_path: str) -> float:
        pass
    
    def get_parent(self, file_path: str) -> str:
        return "dummy"

    def join_path(self, root_path: str, target_path: str) -> str:
        return "dummy"


class MockLoader(AbstractLoader):
    """Mock loader for testing."""

    def __init__(self, filestore: FileStoreInterface, scan_type: RawScanType, score: float):
        """Initialize mock loader with configurable score and scan type."""
        super().__init__(filestore)
        self._scan_type = scan_type
        self._score = score

    def load(self, path: str):
        """Not implemented for testing."""
        pass

    def get_scan_type(self) -> RawScanType:
        """Get the configured scan type."""
        return self._scan_type

    def get_score(self, path: str) -> float:
        """Return the configured score."""
        return self._score

    def parse_metadata(self, path: str):
        """Not implemented for testing."""
        pass
    
    def parse_tavi_metadata(self, path: str) -> dict:
        """Parse metadata."""
        return {}

    def parse_scan_values(self, path: str):
        """Not implemented for testing."""
        pass

    def parse_external_metadata(self, path: str):
        """Not implemented for testing."""
        pass

    def adapt_scan_data(self, meta, values):
        """Not implemented for testing."""
        pass


class TestRawScanClassifier(unittest.TestCase):
    """Test cases for RawScanClassifier."""

    @classmethod
    def setUpClass(cls) -> None:
        """Set up class-level fixtures."""
        # Initialize LoaderRegistry singleton with mock filestore
        cls.mock_filestore = MockFileStore()
        LoaderRegistry().set_filestore(cls.mock_filestore)

    def setUp(self) -> None:
        """Set up test fixtures."""
        self.classifier = RawScanClassifier()

    def test_classifier_initialization(self) -> None:
        """Test that classifier initializes with a loader registry."""
        self.assertIsNotNone(
            self.classifier.loader_registry,
            "Classifier should have a loader registry",
        )

    def test_single_loader_scores_file(self) -> None:
        """Test classification with a single scoring loader."""
        # Mock the registry to return one loader
        mock_loader = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.9,
        )
        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[mock_loader]
        )

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.ORNLSpice,
            "Should return the scan type of the only loader",
        )

    def test_multiple_loaders_highest_score_wins(self) -> None:
        """Test that loader with highest score is selected."""
        # Create multiple loaders with different scores
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.6,
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.9,  # Highest score
        )
        loader3 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.0,
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[loader1, loader2, loader3]
        )

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.ORNLSpice,
            "Should select loader with highest score",
        )

    def test_no_loaders_above_zero_returns_none(self) -> None:
        """Test that NONE is returned when all loaders score 0."""
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.0,
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.0,
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[loader1, loader2]
        )

        result = self.classifier.get_classification("unknown_file")

        self.assertEqual(
            result,
            RawScanType.NONE,
            "Should return NONE when no loader scores above 0",
        )

    def test_first_winner_selected_on_tie(self) -> None:
        """Test that first loader is selected when multiple loaders tie."""
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.8,
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.8,  # Same score as loader1
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[loader1, loader2]
        )

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.ORNLSpice,
            "Should select first loader in case of tie",
        )

    def test_perfect_score_selected_immediately(self) -> None:
        """Test that perfect score (1.0) is recognized."""
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            1.0,  # Perfect score
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.99,
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[loader1, loader2]
        )

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.ORNLSpice,
            "Should select loader with perfect score",
        )

    def test_low_scores_not_selected_over_zero(self) -> None:
        """Test that even very low positive scores beat 0."""
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.01,  # Very low but positive
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.0,
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[loader1, loader2]
        )

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.ORNLSpice,
            "Should select loader with score > 0 over 0",
        )

    def test_classifier_called_with_file_path(self) -> None:
        """Test that loader is called with the correct file path."""
        mock_loader = mock.Mock(spec=AbstractLoader)
        mock_loader.get_score.return_value = 0.5
        mock_loader.get_scan_type.return_value = RawScanType.ORNLSpice

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=[mock_loader]
        )

        test_path = "/path/to/test.dat"
        self.classifier.get_classification(test_path)

        mock_loader.get_score.assert_called_once_with(test_path)

    def test_all_loaders_scored_before_decision(self) -> None:
        """Test that all loaders are scored before selecting winner."""
        loader1 = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.5,
        )
        loader2 = MockLoader(
            self.mock_filestore,
            RawScanType.NONE,
            0.7,
        )

        mock_loaders = [
            mock.Mock(spec=AbstractLoader, get_score=mock.Mock(return_value=0.5),
                     get_scan_type=mock.Mock(return_value=RawScanType.ORNLSpice)),
            mock.Mock(spec=AbstractLoader, get_score=mock.Mock(return_value=0.7),
                     get_scan_type=mock.Mock(return_value=RawScanType.NONE)),
        ]

        self.classifier.loader_registry.get_loaders = mock.Mock(
            return_value=mock_loaders
        )

        self.classifier.get_classification("test.dat")

        # Verify both loaders were scored
        for mock_loader in mock_loaders:
            mock_loader.get_score.assert_called_once()

    def test_score_range_zero_to_one(self) -> None:
        """Test that scores between 0.0 and 1.0 are handled correctly."""
        scores = [0.0, 0.25, 0.5, 0.75, 1.0]
        scan_types = [
            RawScanType.NONE,
            RawScanType.ORNLSpice,
            RawScanType.NONE,
            RawScanType.ORNLSpice,
            RawScanType.NONE,
        ]

        for score, scan_type in zip(scores, scan_types):
            loader = MockLoader(
                self.mock_filestore,
                scan_type,
                score,
            )
            self.classifier.loader_registry.get_loaders = mock.Mock(
                return_value=[loader]
            )

            result = self.classifier.get_classification("test.dat")

            if score == 0.0:
                self.assertEqual(
                    result,
                    RawScanType.NONE,
                    f"Score {score} should result in NONE",
                )
            else:
                self.assertEqual(
                    result,
                    scan_type,
                    f"Score {score} should be accepted",
                )

    def test_empty_loader_list_returns_none(self) -> None:
        """Test that NONE is returned when no loaders are registered."""
        self.classifier.loader_registry.get_loaders = mock.Mock(return_value=[])

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.NONE,
            "Should return NONE when no loaders are registered",
        )

    def test_initial_top_pick_is_none_zero(self) -> None:
        """Test that initial tracking starts with NONE and 0 score."""
        loader = MockLoader(
            self.mock_filestore,
            RawScanType.ORNLSpice,
            0.0,  # Score below initial
        )

        self.classifier.loader_registry.get_loaders = mock.Mock(return_value=[loader])

        result = self.classifier.get_classification("test.dat")

        self.assertEqual(
            result,
            RawScanType.NONE,
            "Score of 0 should not beat initial NONE",
        )


if __name__ == "__main__":
    unittest.main()
