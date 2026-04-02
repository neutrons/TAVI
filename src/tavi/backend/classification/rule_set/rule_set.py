"""Rule set for organizing classification rules."""

import math

from tavi.backend.classification.rule.interface.rule_interface import RuleInterface


class RuleSet:
    """A collection of classification rules with weights."""

    def __init__(self) -> None:
        """Initialize an empty rule set."""
        self.rules: dict[RuleInterface, int] = {}

    def register(self, rule: RuleInterface, weight: int) -> None:
        """
        Register a rule with a weight.

        Args:
            rule: The rule to register.
            weight: The weight of the rule.

        """
        self.rules[rule] = weight

    def validate(self) -> None:
        """Validate that weights total to 1.0."""
        total_weight = 0
        for rule in self.get_rules():
            total_weight += self.get_weight(rule)
        if not math.isclose(total_weight, 1):
            raise ValueError("RuleSet incorrectly configured.  Weights should total to ~1")

    def get_rules(self) -> list[RuleInterface]:
        """
        Get all registered rules.

        Returns:
            A list of registered rules.

        """
        return list(self.rules.keys())

    def get_weight(self, rule: RuleInterface) -> int:
        """
        Get the weight of a rule.

        Args:
            rule: The rule to get the weight for.

        Returns:
            The weight of the rule.

        """
        return self.rules[rule]
