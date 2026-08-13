RuleBasedClassifier
===================

Overview
--------

The ``RuleBasedClassifier`` is a flexible classification system that uses a set of rules to score and classify files. It's designed to determine file types or characteristics by applying multiple scoring rules and combining their results.

Key Features:

- **Rule-based scoring**: Assigns scores to files based on configurable rules
- **Weighted rules**: Each rule can have a different weight to influence the final score
- **Extensible**: Easy to create custom rules by implementing the ``RuleInterface``
- **File store integration**: Uses a file store interface to access file properties and contents

Architecture
------------

The RuleBasedClassifier is built on three core components:

**RuleBasedClassifier**
    The main classifier that applies a set of rules to score files. It computes a
    cumulative score by multiplying each rule's score by its weight. It is a
    ``@Singleton``, so the ``filestore`` constructor argument is only honoured on
    the first construction; later callers rebind it with ``set_filestore()``
    (which is what ``ORNLSpiceLoader.get_score`` does on every call).

**RuleInterface**
    An abstract base class for all classification rules. Each rule implements a single classification criterion.

    .. code-block:: python

        def get_score(self, path: str, filestore: FileStoreInterface) -> int:
            """Get the score for a file path."""

    The interface declares ``int``, and every built-in rule returns ``0`` or
    ``1``. Custom rules may return a fractional confidence — the classifier sums
    into a float — but the declared type has not been widened, so treat 0/1 as
    the convention.

**RuleSet**
    A collection of rules with associated weights. Each rule is registered with a weight that determines its influence on the final score.

How It Works
~~~~~~~~~~~~

1. Create or retrieve a ``RuleSet`` containing registered rules and their weights
2. Pass the file path and ``RuleSet`` to the classifier's ``get_score()`` method
3. The classifier iterates through all rules in the set
4. For each rule, it calculates: ``rule_score × rule_weight``
5. Returns the cumulative sum of all weighted rule scores

The higher the final score, the better the file matches the criteria defined by the rules.

Basic Usage
-----------

.. code-block:: python

    from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
    from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
    from tavi.library.storage.local_file_store import LocalFileStore

    # Initialize the file store
    filestore = LocalFileStore()

    # Create a classifier (a singleton - the filestore is only bound on first call)
    classifier = RuleBasedClassifier(filestore)

    # Create a rule set
    rule_set = ORNLSpiceRuleSet()

    # Score a file
    score = classifier.get_score("/path/to/file.dat", rule_set)
    print(f"Classification score: {score}")

Creating Custom Rules
---------------------

To create a custom rule, inherit from ``RuleInterface`` and implement the ``get_score()`` method:

.. code-block:: python

    from tavi.backend.classification.rule.interface.rule_interface import RuleInterface
    from tavi.library.storage.interface.file_store_interface import FileStoreInterface

    class MyCustomRule(RuleInterface):
        """Custom rule that checks for specific file characteristics."""

        def get_score(self, path: str, filestore: FileStoreInterface) -> float:
            """
            Score a file based on custom criteria.

            Args:
                path: The file path to evaluate.
                filestore: The file store interface to access file properties.

            Returns:
                The score as a float between 0.0 and 1.0.

            """
            # Example: Check if file has a specific extension
            if path.endswith('.txt'):
                return 1.0
            return 0.0

Creating Custom Rule Sets
--------------------------

To create a custom rule set, inherit from ``RuleSet`` and register your rules. **Important: All weights must sum to 1.0.**

.. code-block:: python

    from tavi.backend.classification.rule_set.rule_set import RuleSet
    from your_module import MyCustomRule, AnotherRule

    class MyRuleSet(RuleSet):
        """Custom rule set for specialized file classification."""

        def __init__(self) -> None:
            """Initialize the custom rule set with registered rules."""
            super().__init__()
            # Register rules with weights that total to 1.0
            self.register(MyCustomRule(), 0.6)      # Weight of 0.6
            self.register(AnotherRule(), 0.4)       # Weight of 0.4
            # Validate that weights sum to 1.0
            self.validate()

Built-in Rule Sets
-------------------

**ORNLSpiceRuleSet**
    A pre-configured rule set for classifying ORNL SPICE format files. Five rules,
    each registered with weight ``1/5`` so the total is 1.0 and a perfectly
    matching file scores 1.0:

    ``InstrumentInFilenameRule``
        The filename's leading ``_``-delimited token, upper-cased, is in
        ``Config["ORNL.instrument.names"]`` (``HB1A``, ``HB1``, ``HB3``, ``CG4C``).

    ``DefXYRule``
        Every field in ``Config["ORNL.spice.metadata.field-name.def"]``
        (``def_x``, ``def_y``) appears as a header line of the form
        ``# <field> =`` followed by a space.

    ``DatFormatRule``
        The file parses as whitespace-separated columns with at least one numeric
        column, either treating ``#`` as a comment marker or starting from the
        first all-numeric line.

    ``HashtagCommentRule``
        At least one line begins with ``#``.

    ``SpiceFileExtRule``
        The file extension equals ``Config["ORNL.spice.file.ext"]`` (``.dat``).

    Every rule swallows its own exceptions and returns ``0``, so an unreadable or
    malformed file scores low rather than raising during classification.

Example: ORNL SPICE Classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.backend.classification.rule_based_classifier import RuleBasedClassifier
    from tavi.backend.classification.rule_set.ornl_spice_rule_set import ORNLSpiceRuleSet
    from tavi.library.storage.local_file_store import LocalFileStore

    filestore = LocalFileStore()
    classifier = RuleBasedClassifier(filestore)
    rule_set = ORNLSpiceRuleSet()

    # Score multiple files
    files = [
        "HB1A_exp0978_scan0001.dat",
        "sample.txt",
        "data.spice",
    ]

    for filepath in files:
        score = classifier.get_score(filepath, rule_set)
        if score > 0.5:  # Example threshold
            print(f"{filepath} is likely an ORNL SPICE file (score: {score:.2f})")

Implementation Details
----------------------

**Score Aggregation**
    The classifier accumulates scores from all rules:

    .. code-block:: python

        score = 0
        for rule in rule_set.get_rules():
            score += rule.get_score(path, filestore) * rule_set.get_weight(rule)

**File Store Integration**
    Rules access file properties through the ``FileStoreInterface``, which allows rules to:

    - Read file contents (``read_text_file``)
    - Check file metadata — ``get_file_name``, ``get_file_ext``,
      ``get_file_size_mb``
    - Validate file format (``validate_file``)

**Extensibility**
    The design supports:

    - Custom rules that implement any classification logic
    - Custom rule sets with different weight configurations
    - Weighted aggregation for fine-tuned classification

Best Practices
--------------

1. **Rule Design**: Keep rules simple and focused on a single classification criterion
2. **Weight Configuration**: Use weights to prioritize important rules over less critical ones
3. **Score Interpretation**: Define thresholds for your use case (e.g., score > 0.5 means "likely a match"). Since scores are between 0.0 and 1.0, consider what confidence level you require.
4. **Testing**: Test your custom rules with representative file samples
5. **Documentation**: Document the purpose and expected behavior of custom rules
