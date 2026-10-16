"""Filter view widget."""

from collections.abc import Callable
from typing import Any
from uuid import uuid4

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QComboBox,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)


class CriteriaItemWidget(QWidget):
    """
    Criteria Item Widget.

    Holds atomic filtering statements.
    """

    operators = ["is", "is not", "has", "has not", "=", "!=", ">", "<", ">=", "<="]

    def __init__(self, parent: QWidget | None = None) -> None:
        """Construct individual criteria."""
        super().__init__(parent)
        criteria_layout = QHBoxLayout(self)

        self.uuid = uuid4()
        self.criteria_key = QLineEdit("key")
        self.criteria_operator = QComboBox(editable=False)
        self.criteria_operator.addItems(self.operators)
        self.criteria_value = QLineEdit("value")
        self.criteria_delete_button = QPushButton()

        style = self.style()
        trash_icon = style.standardIcon(style.StandardPixmap.SP_TrashIcon)
        self.criteria_delete_button.setIcon(trash_icon)

        criteria_layout.addWidget(self.criteria_key)
        criteria_layout.addWidget(self.criteria_operator)
        criteria_layout.addWidget(self.criteria_value)
        criteria_layout.addWidget(self.criteria_delete_button)

    def set_sample(self) -> None:
        """Set criteria to sample mode: basic values and read-only."""
        self.criteria_key.setText("scan")
        self.criteria_key.setEnabled(False)
        self.criteria_operator.clear()
        self.criteria_operator.addItem(">=")
        self.criteria_operator.setEnabled(False)
        self.criteria_value.setText("0")
        self.criteria_value.setEnabled(False)
        self.criteria_delete_button.setEnabled(False)


class CriteriaSetScrollArea(QScrollArea):
    """
    Criteria Set Widget.

    Holds a set of filter statements to create a filter.
    Handles adding and removing criteria from the set.
    """

    def __init__(self, parent: QWidget | None = None, set_sample: bool = False) -> None:
        """Construct criteria set object."""
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        dummy_widget = QWidget()
        self.setWidget(dummy_widget)
        self.criteria_list_layout = QVBoxLayout(dummy_widget)
        self.criteria_list_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.criteria_store = dict()
        initial_criteria = CriteriaItemWidget()

        self.criteria_store[initial_criteria.uuid.hex] = initial_criteria
        self.criteria_list_layout.addWidget(initial_criteria)
        self.add_criteria_button = QPushButton("+")
        if not set_sample:
            initial_criteria.criteria_delete_button.clicked.connect(self.delete_criteria(initial_criteria))

            self.add_criteria_button.clicked.connect(self.add_criteria)

        self.criteria_list_layout.addWidget(self.add_criteria_button)
        if set_sample:
            initial_criteria.set_sample()
            self.add_criteria_button.setEnabled(False)

    def add_criteria(self) -> None:
        """Add an empty criteria to the set."""
        new_criteria = CriteriaItemWidget()
        self.criteria_store[new_criteria.uuid.hex] = new_criteria
        self.criteria_list_layout.addWidget(new_criteria)
        new_criteria.criteria_delete_button.clicked.connect(self.delete_criteria(new_criteria))
        self.criteria_list_layout.addWidget(self.add_criteria_button)
        # Re-adding an existing widget just moves it to the end

    def delete_criteria(self, criteria: CriteriaItemWidget) -> Callable:
        """Create a callable to delete a provided criteria."""

        def _delete_criteria() -> None:
            """Hides and deletes the externally provided criteria from the set."""
            criteria.hide()
            self.criteria_store.pop(criteria.uuid.hex)
            self.criteria_list_layout.removeWidget(criteria)
            criteria.deleteLater()

        return _delete_criteria


class FilterView(QWidget):
    """Filter panel widget."""

    def __init__(self, parent: Any = None) -> None:
        """Construct filter view."""
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the filter UI."""
        layout = QVBoxLayout(self)
        filter_box = QGroupBox()
        # Box for a nice outline
        filter_components_layout = QVBoxLayout()

        # ---
        # Upper Filter Controls
        # ---
        filter_header_layout = QHBoxLayout()

        filter_applied_layout = QHBoxLayout()
        filter_applied_layout.setDirection(filter_applied_layout.Direction.LeftToRight)
        # Default but provided for posterity

        # Start with the indicator for current filter applied
        filter_applied_prefix = QLabel("Current Filter:")
        filter_applied_prefix.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)
        filter_applied_layout.addWidget(filter_applied_prefix)

        self.filter_applied_friendly_name = QLabel(" " * 5)
        # This will react to the currently applied filter here or elsewhere
        # self.filter_applied_friendly_name.setWordWrap(True)
        self.filter_applied_friendly_name.setFrameStyle(QFrame.Shape.Box)
        filter_applied_layout.addWidget(self.filter_applied_friendly_name)

        filter_header_layout.addLayout(filter_applied_layout)

        # Next the new, del, and dropdown
        filter_controls_layout = QHBoxLayout()
        filter_controls_layout.setDirection(filter_controls_layout.Direction.RightToLeft)

        self.filter_dropdown = QComboBox(editable=True, insertPolicy=QComboBox.InsertPolicy.InsertAtCurrent)
        self.filter_dropdown.addItem("Sample", uuid4())
        self.filter_dropdown.editTextChanged.connect(self.immutable_sample_name)
        filter_controls_layout.addWidget(self.filter_dropdown)

        self.filter_del_button = QPushButton("Del")
        self.filter_del_button.clicked.connect(self.del_filter)
        self.filter_del_button.setEnabled(False)
        filter_controls_layout.addWidget(self.filter_del_button)

        self.filter_new_button = QPushButton("New")
        self.filter_new_button.pressed.connect(self.new_filter)
        filter_controls_layout.addWidget(self.filter_new_button)

        filter_controls_label = QLabel("Filters")
        filter_controls_layout.addWidget(filter_controls_label)

        filter_header_layout.addLayout(filter_controls_layout)

        # ---
        # Criteria Controls
        # ---
        self.filter_criteria_set_list = QStackedWidget()
        self.filter_criteria_set_list.addWidget(CriteriaSetScrollArea(set_sample=True))
        self.filter_dropdown.currentIndexChanged.connect(self.propagate_index_change)

        # ---
        # Filter Application Controls
        # ---
        self.applied_filter = None
        filter_footer_layout = QHBoxLayout()

        filter_left_footer = QHBoxLayout()
        self.clear_filter_button = QPushButton("Clear Filter")
        self.clear_filter_button.clicked.connect(self.clear_filter)
        filter_left_footer.addWidget(self.clear_filter_button)
        filter_footer_layout.addLayout(filter_left_footer)

        filter_right_footer = QHBoxLayout()
        filter_right_footer.setDirection(filter_right_footer.Direction.RightToLeft)
        self.apply_filter_button = QPushButton("Apply Filter")
        self.apply_filter_button.clicked.connect(self.apply_filter)
        filter_right_footer.addWidget(self.apply_filter_button)
        filter_footer_layout.addLayout(filter_right_footer)

        # ---
        # Combine Main Components
        # ---
        filter_components_layout.addLayout(filter_header_layout)
        filter_components_layout.addWidget(self.filter_criteria_set_list)
        filter_components_layout.addLayout(filter_footer_layout)
        filter_box.setLayout(filter_components_layout)
        layout.addWidget(filter_box)

    def new_filter(self) -> None:
        """Create a new filter in the view."""
        filter_index = self.filter_dropdown.count()
        self.filter_dropdown.addItem(f"FilterName{filter_index}", uuid4())
        self.filter_criteria_set_list.addWidget(CriteriaSetScrollArea())
        self.filter_dropdown.setCurrentIndex(self.filter_dropdown.count() - 1)

    def del_filter(self) -> None:
        """Delete the currently highlighted filter."""
        index_for_deletion = self.filter_dropdown.currentIndex()
        if index_for_deletion == 0:
            return
        if self.filter_dropdown.itemData(index_for_deletion) == self.applied_filter:
            self.clear_filter()
        self.filter_dropdown.removeItem(index_for_deletion)
        set_for_deletion = self.filter_criteria_set_list.widget(index_for_deletion)
        if set_for_deletion is not None:
            self.filter_criteria_set_list.removeWidget(set_for_deletion)
        self.filter_dropdown.setCurrentIndex(0)

    def immutable_sample_name(self) -> None:
        """Revert any attempt to modify the sample filter name."""
        sample_text = self.filter_dropdown.itemText(0)
        if sample_text != "Sample":
            self.filter_dropdown.setItemText(0, "Sample")

    def propagate_index_change(self) -> None:
        """Propagate relevant changes from selecting a different filter index."""
        dropdown_index = self.filter_dropdown.currentIndex()
        self.filter_del_button.setEnabled(True)
        if dropdown_index == 0:
            self.filter_del_button.setEnabled(False)
        self.filter_criteria_set_list.setCurrentIndex(dropdown_index)

    def apply_filter(self) -> None:
        """Apply the currently highlighted filter with appropriate messaging."""
        # Do bundling and sending of filter things here
        self.applied_filter = self.filter_dropdown.currentData()
        # Replace below with basing the label on received filter event
        friendly_name = self.filter_dropdown.currentText()
        if len(friendly_name) > 30:
            friendly_name = friendly_name[:30] + "..."
        self.filter_applied_friendly_name.setText(friendly_name)

    def clear_filter(self) -> None:
        """Clear the currently applied filter with appropriate messaging."""
        # Do sending of clearing filter here
        self.applied_filter = None
        self.filter_applied_friendly_name.setText(" " * 5)
