"""Raw scan classifier for identifying scan file types."""

import logging

from tavi.library.data.enum.raw_scan_type import RawScanType
from tavi.library.storage.loader.interface.base import AbstractLoader
from tavi.library.storage.loader.loader_registry import LoaderRegistry

logger = logging.getLogger(__name__)


class RawScanClassifier:
    """Classifies raw scan files using registered loaders."""

    def __init__(self) -> None:
        """Initialize classifier with loader registry."""
        self.loader_registry: LoaderRegistry = LoaderRegistry()

    def get_classification(self, file_path: str) -> RawScanType:
        """Get classification for a file by score from all loaders."""
        loaders: list[AbstractLoader] = self.loader_registry.get_loaders()
        top_pick: tuple[RawScanType, float] = (RawScanType.NONE, 0)
        for loader in loaders:
            score = loader.get_score(file_path)
            if top_pick[1] < score:
                top_pick = (loader.get_scan_type(), score)
            logger.debug(f"Loader {loader.get_scan_type()} rated file {file_path} with score: {score}.")
        return top_pick[0]
