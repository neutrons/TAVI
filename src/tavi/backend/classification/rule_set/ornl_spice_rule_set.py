"""ORNL SPICE file rule set."""

from tavi.backend.classification.rule.dat_format_rule import DatFormatRule
from tavi.backend.classification.rule.hashtag_comment_rule import HashtagCommentRule
from tavi.backend.classification.rule.ORNL.def_xy_rule import DefXYRule
from tavi.backend.classification.rule.ORNL.inst_in_filename_rule import (
    InstrumentInFilenameRule,
)
from tavi.backend.classification.rule.ORNL.spice_file_ext_rule import SpiceFileExtRule
from tavi.backend.classification.rule_set.rule_set import RuleSet


class ORNLSpiceRuleSet(RuleSet):
    """Rule set for ORNL SPICE files."""

    def __init__(self) -> None:
        """Initialize the ORNL SPICE rule set with default rules."""
        super().__init__()
        self.register(InstrumentInFilenameRule(), 1)
        self.register(DefXYRule(), 1)
        self.register(DatFormatRule(), 1)
        self.register(HashtagCommentRule(), 1)
        self.register(SpiceFileExtRule(), 1)
