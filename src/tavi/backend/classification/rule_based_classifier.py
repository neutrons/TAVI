"""Rule-based classifier for files."""

from tavi.backend.classification.rule_set.rule_set import RuleSet
from tavi.library.storage.interface.file_store_interface import FileStoreInterface


class RuleBasedClassifier:
    """Classifier that uses a set of rules to score files."""

    def __init__(self, filestore: FileStoreInterface) -> None:
        """
        Initialize the classifier with a file store.

        Args:
            filestore: The file store interface to use.

        """
        self.filestore: FileStoreInterface = filestore

    def _get_score_from_rules(self, path: str, rule_set: RuleSet) -> float:
        """
        Get the cumulative score from all rules.

        Args:
            path: The file path to score.
            rule_set: The set of rules to apply.

        Returns:
            The cumulative score.

        """
        score: float = 0
        for rule in rule_set.get_rules():
            score += rule.get_score(path, self.filestore) * rule_set.get_weight(rule)
        return score

    def get_score(self, path: str, rule_set: RuleSet) -> int:
        """
        Get the score for a file path.

        Args:
            path: The file path to score.
            rule_set: The set of rules to apply.

        Returns:
            The score for the file path.

        """
        return self._get_score_from_rules(path, rule_set)
