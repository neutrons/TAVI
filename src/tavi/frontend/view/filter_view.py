"""Filter view widget."""

from collections.abc import Callable
from typing import Any

from neutrons_standard.config import Config
from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QComboBox,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from tavi.library.data.scan import UUID, new_uuid

Criterion = tuple[str, str, str]
"""One criterion as the panel shows it: (key, operator, value)."""


class CriteriaItemWidget(QWidget):
    """
    Criteria Item Widget.

    Holds atomic filtering statements.
    """

    changed = Signal()
    operators = ["is", "is not", "has", "has not", "=", "!=", ">", "<", ">=", "<="]

    def __init__(self, parent: QWidget | None = None) -> None:
        """Construct individual criteria."""
        super().__init__(parent)
        criteria_layout = QHBoxLayout(self)

        self.uuid = new_uuid()
        self.criteria_key = QLineEdit()
        self.criteria_key.setPlaceholderText("key")
        self.criteria_operator = QComboBox(editable=False)
        self.criteria_operator.addItems(self.operators)
        self.criteria_value = QLineEdit()
        self.criteria_value.setPlaceholderText("value")
        self.criteria_delete_button = QPushButton()

        style = self.style()
        trash_icon = style.standardIcon(style.StandardPixmap.SP_TrashIcon)
        self.criteria_delete_button.setIcon(trash_icon)

        criteria_layout.addWidget(self.criteria_key)
        criteria_layout.addWidget(self.criteria_operator)
        criteria_layout.addWidget(self.criteria_value)
        criteria_layout.addWidget(self.criteria_delete_button)

        self.criteria_key.textChanged.connect(self.changed)
        self.criteria_operator.currentTextChanged.connect(self.changed)
        self.criteria_value.textChanged.connect(self.changed)

    def set_criterion(self, criterion: Criterion) -> None:
        """Show ``criterion`` in the key, operator and value fields."""
        key, operator, value = criterion
        self.criteria_key.setText(key)
        self.criteria_operator.setCurrentText(operator)
        self.criteria_value.setText(value)

    def criterion(self) -> Criterion:
        """Return the key, operator and value as currently entered."""
        return self.criteria_key.text(), self.criteria_operator.currentText(), self.criteria_value.text()

    def set_sample(self) -> None:
        """Set criteria to sample mode: basic values and read-only."""
        self.criteria_key.setText(Config["filter.sample.criterion.key"])
        self.criteria_key.setEnabled(False)
        self.criteria_operator.clear()
        self.criteria_operator.addItem(Config["filter.sample.criterion.operator"])
        self.criteria_operator.setEnabled(False)
        self.criteria_value.setText(str(Config["filter.sample.criterion.value"]))
        self.criteria_value.setEnabled(False)
        self.criteria_delete_button.setEnabled(False)


class CriteriaSetScrollArea(QScrollArea):
    """
    Criteria Set Widget.

    Holds a set of filter statements to create a filter, and the filter's name as currently entered.
    Handles adding and removing criteria from the set, and remembers what the filter was last saved as.
    """

    changed = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        set_sample: bool = False,
        name: str = "",
        criteria: list[Criterion] | None = None,
        saved: bool = True,
    ) -> None:
        """
        Construct a criteria set showing ``criteria``, or one blank criterion if none are given.

        ``saved`` says whether the filter is in the project as shown; if not, it stays unsaved until ``mark_saved``.
        """
        super().__init__(parent)
        self.name = Config["filter.sample.name"] if set_sample else name
        self._saved: tuple[str, list[Criterion]] | None = None
        self.setWidgetResizable(True)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        dummy_widget = QWidget()
        self.setWidget(dummy_widget)
        self.criteria_list_layout = QVBoxLayout(dummy_widget)
        self.criteria_list_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.criteria_store = dict()
        self.add_criteria_button = QPushButton("+")
        self.criteria_list_layout.addWidget(self.add_criteria_button)

        if set_sample:
            sample_criteria = CriteriaItemWidget()
            sample_criteria.set_sample()
            self.criteria_store[sample_criteria.uuid] = sample_criteria
            self.criteria_list_layout.insertWidget(0, sample_criteria)
            self.add_criteria_button.setEnabled(False)
            self.mark_saved(self.name, self.criteria())
            return

        self.add_criteria_button.clicked.connect(self.add_criteria)
        for criterion in criteria if criteria is not None else [("", "is", "")]:
            self.add_criteria().set_criterion(criterion)
        if saved:
            self.mark_saved(self.name, self.criteria())

    def add_criteria(self) -> CriteriaItemWidget:
        """Add an empty criteria to the set and return it."""
        new_criteria = CriteriaItemWidget()
        self.criteria_store[new_criteria.uuid] = new_criteria
        self.criteria_list_layout.addWidget(new_criteria)
        new_criteria.criteria_delete_button.clicked.connect(self.delete_criteria(new_criteria))
        new_criteria.changed.connect(self.changed)
        self.criteria_list_layout.addWidget(self.add_criteria_button)
        # Re-adding an existing widget just moves it to the end
        self.changed.emit()
        return new_criteria

    def criteria(self) -> list[Criterion]:
        """Return every criterion in the set, in the order shown."""
        return [item.criterion() for item in self.criteria_store.values()]

    def mark_saved(self, name: str, criteria: list[Criterion]) -> None:
        """Remember that the filter was saved with ``name`` and ``criteria``."""
        self._saved = (name, [tuple(criterion) for criterion in criteria])
        self.changed.emit()

    def is_in_project(self) -> bool:
        """Return whether the filter has been saved into the project at all."""
        return self._saved is not None

    def is_unsaved(self) -> bool:
        """Return whether the filter was never saved, or its name or criteria differ from what it was saved as."""
        return self._saved != (self.name, self.criteria())

    def delete_criteria(self, criteria: CriteriaItemWidget) -> Callable:
        """Create a callable to delete a provided criteria."""

        def _delete_criteria() -> None:
            """Hides and deletes the externally provided criteria from the set."""
            criteria.hide()
            self.criteria_store.pop(criteria.uuid)
            self.criteria_list_layout.removeWidget(criteria)
            criteria.deleteLater()
            self.changed.emit()

        return _delete_criteria


class FilterView(QWidget):
    """
    Filter panel widget.

    Shows the filters ``FilterPresenter`` tells it about; its buttons only report what the user asked for.
    Emit the ``*_signal`` signals rather than calling their slots, so updates published from a model's
    worker thread land on the GUI thread.
    """

    DEFAULT_APPLIED_TEXT = "Default Applied"
    UNSAVED_MARKER = "*"

    add_filter_signal = Signal(object, str, list)
    sync_filter_signal = Signal(object, str, list)
    remove_filter_signal = Signal(object)
    show_applied_filter_signal = Signal(str)
    show_default_applied_signal = Signal()
    show_validation_errors_signal = Signal(list)

    def __init__(self, parent: Any = None) -> None:
        """Construct filter view."""
        super().__init__(parent)
        self._build_ui()
        # AutoConnection: direct call on the GUI thread (tests), queued hop when emitted from
        # a worker thread (FilterModel/TaviProjectModel running behind their Proxy).
        self.add_filter_signal.connect(self._add_filter)
        self.sync_filter_signal.connect(self._sync_filter)
        self.remove_filter_signal.connect(self._remove_filter)
        self.show_applied_filter_signal.connect(self._show_applied_filter)
        self.show_default_applied_signal.connect(self._show_default_applied)
        self.show_validation_errors_signal.connect(self._show_validation_errors)

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

        self.filter_applied_friendly_name = QLabel(self.DEFAULT_APPLIED_TEXT)
        # This will react to the currently applied filter here or elsewhere
        # self.filter_applied_friendly_name.setWordWrap(True)
        self.filter_applied_friendly_name.setFrameStyle(QFrame.Shape.Box)
        # Greyed out while no filter is applied, so the default doesn't read as a saved filter's name.
        self.filter_applied_friendly_name.setEnabled(False)
        filter_applied_layout.addWidget(self.filter_applied_friendly_name)

        filter_header_layout.addLayout(filter_applied_layout)

        # Next the new, del, and dropdown
        filter_controls_layout = QHBoxLayout()
        filter_controls_layout.setDirection(filter_controls_layout.Direction.RightToLeft)

        # The line edit renames the picked filter as the user types (see _rename_current_filter), so Enter
        # needn't insert anything, and must not jump to another filter that happens to share the name.
        self.filter_dropdown = QComboBox(editable=True, insertPolicy=QComboBox.InsertPolicy.NoInsert)
        self.filter_dropdown.setDuplicatesEnabled(True)
        self.filter_dropdown.setCompleter(None)
        self.filter_dropdown.addItem(Config["filter.sample.name"], UUID(value=Config["filter.sample.uuid"]))
        self.filter_dropdown.editTextChanged.connect(self.immutable_sample_name)
        self.filter_dropdown.lineEdit().textEdited.connect(self._rename_current_filter)
        filter_controls_layout.addWidget(self.filter_dropdown)

        self.filter_del_button = QPushButton("Del")
        self.filter_del_button.setEnabled(False)
        filter_controls_layout.addWidget(self.filter_del_button)

        self.filter_new_button = QPushButton("New")
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
        filter_footer_layout = QHBoxLayout()

        filter_left_footer = QHBoxLayout()
        self.clear_filter_button = QPushButton("Clear Filter")
        filter_left_footer.addWidget(self.clear_filter_button)
        filter_footer_layout.addLayout(filter_left_footer)

        filter_right_footer = QHBoxLayout()
        filter_right_footer.setDirection(filter_right_footer.Direction.RightToLeft)
        self.apply_filter_button = QPushButton("Save && Apply Filter")
        filter_right_footer.addWidget(self.apply_filter_button)
        filter_footer_layout.addLayout(filter_right_footer)

        # ---
        # Combine Main Components
        # ---
        filter_components_layout.addLayout(filter_header_layout)
        filter_components_layout.addWidget(self.filter_criteria_set_list)
        self.filter_validation_errors = QLabel()
        self.filter_validation_errors.setWordWrap(True)
        self.filter_validation_errors.setStyleSheet("color: red")
        self.filter_validation_errors.hide()
        filter_components_layout.addWidget(self.filter_validation_errors)
        filter_components_layout.addLayout(filter_footer_layout)
        filter_box.setLayout(filter_components_layout)
        layout.addWidget(filter_box)

    def current_filter_uuid(self) -> UUID:
        """Return the uuid of the filter picked in the dropdown."""
        return self.filter_dropdown.currentData()

    def current_filter_name(self) -> str:
        """Return the picked filter's name as currently entered, without the unsaved-changes marker."""
        return self.filter_criteria_set_list.currentWidget().name

    def current_criteria(self) -> list[Criterion]:
        """Return the picked filter's criteria as currently entered."""
        return self.filter_criteria_set_list.currentWidget().criteria()

    def confirm_delete_filter(self) -> bool:
        """Ask the user to confirm deleting the picked filter; return whether they did."""
        reply = QMessageBox.question(
            self,
            "Delete Filter",
            f"Delete filter '{self.current_filter_name()}'? This cannot be undone.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return reply == QMessageBox.StandardButton.Yes

    def filter_count(self) -> int:
        """Return how many filters the dropdown lists, Sample included."""
        return self.filter_dropdown.count()

    def current_filter_in_project(self) -> bool:
        """Return whether the filter picked in the dropdown has been saved into the project."""
        return self.filter_criteria_set_list.widget(self.filter_dropdown.currentIndex()).is_in_project()

    def add_unsaved_filter(self, uuid: UUID, name: str, criteria: list[Criterion]) -> None:
        """List a filter that is not in the project yet, marked unsaved, and pick it."""
        self._list_filter(uuid, name, criteria, saved=False)

    def _add_filter(self, uuid: UUID, name: str, criteria: list[Criterion]) -> None:
        """List a filter newly saved to the project and pick it. Use ``add_filter_signal``."""
        index = self._index_of(uuid)
        if index > 0:
            # Listed already, as an unsaved filter that has now been saved for the first time.
            self.filter_criteria_set_list.widget(index).mark_saved(name, criteria)
            return
        self._list_filter(uuid, name, criteria, saved=True)

    def _list_filter(self, uuid: UUID, name: str, criteria: list[Criterion], saved: bool) -> None:
        """Add a filter to the dropdown and the criteria sets, and pick it."""
        self.filter_dropdown.addItem(name, uuid)
        criteria_set = CriteriaSetScrollArea(name=name, criteria=criteria, saved=saved)
        self.filter_criteria_set_list.addWidget(criteria_set)
        # Listed before the connection, so an unsaved filter gets its marker here.
        self._show_unsaved_marker(criteria_set)
        criteria_set.changed.connect(lambda: self._show_unsaved_marker(criteria_set))
        self.filter_dropdown.setCurrentIndex(self.filter_dropdown.count() - 1)

    def _sync_filter(self, uuid: UUID, name: str, criteria: list[Criterion]) -> None:
        """Remember what a filter was saved as, to compare the panel against. Use ``sync_filter_signal``."""
        index = self._index_of(uuid)
        if index <= 0:
            return
        self.filter_criteria_set_list.widget(index).mark_saved(name, criteria)

    def _index_of(self, uuid: UUID) -> int:
        """Return the dropdown index listing the filter ``uuid``, or -1 if none does."""
        # Not findData: Qt compares wrapped Python objects by identity, and this uuid is an event's copy.
        for index in range(self.filter_dropdown.count()):
            if self.filter_dropdown.itemData(index) == uuid:
                return index
        return -1

    def _rename_current_filter(self, text: str) -> None:
        """Take what the user typed as the picked filter's new, unsaved name."""
        criteria_set = self.filter_criteria_set_list.currentWidget()
        if self.filter_dropdown.currentIndex() == 0:
            self.filter_dropdown.lineEdit().setText(criteria_set.name)
            return
        # The marker shares the line edit with the name, so names can't contain it.
        criteria_set.name = text.replace(self.UNSAVED_MARKER, "")
        self._show_unsaved_marker(criteria_set)

    def _show_unsaved_marker(self, criteria_set: CriteriaSetScrollArea) -> None:
        """List a filter by its name, with the unsaved-changes marker if it differs from what was saved."""
        index = self.filter_criteria_set_list.indexOf(criteria_set)
        if index < 0:
            return
        text = criteria_set.name + (self.UNSAVED_MARKER if criteria_set.is_unsaved() else "")
        if self.filter_dropdown.itemText(index) == text:
            return
        line_edit = self.filter_dropdown.lineEdit()
        cursor = line_edit.cursorPosition()
        # Resets the line edit's text, and its cursor, when this is the picked filter.
        self.filter_dropdown.setItemText(index, text)
        if index == self.filter_dropdown.currentIndex():
            line_edit.setCursorPosition(min(cursor, len(criteria_set.name)))

    def _remove_filter(self, uuid: UUID) -> None:
        """Stop listing a filter, falling back to Sample. Use ``remove_filter_signal``."""
        index = self._index_of(uuid)
        if index <= 0:
            return
        # Pick Sample first, so removing the item doesn't briefly pick whichever filter sits next to it.
        self.filter_dropdown.setCurrentIndex(0)
        self.filter_dropdown.removeItem(index)
        criteria_set = self.filter_criteria_set_list.widget(index)
        if criteria_set is not None:
            self.filter_criteria_set_list.removeWidget(criteria_set)
            criteria_set.deleteLater()

    def immutable_sample_name(self) -> None:
        """Revert any attempt to modify the sample filter name."""
        sample_name = Config["filter.sample.name"]
        if self.filter_dropdown.itemText(0) != sample_name:
            self.filter_dropdown.setItemText(0, sample_name)

    def propagate_index_change(self) -> None:
        """Propagate relevant changes from selecting a different filter index."""
        dropdown_index = self.filter_dropdown.currentIndex()
        self.filter_del_button.setEnabled(True)
        if dropdown_index == 0:
            self.filter_del_button.setEnabled(False)
        self.filter_criteria_set_list.setCurrentIndex(dropdown_index)

    def _show_applied_filter(self, name: str) -> None:
        """Name the applied filter. Use ``show_applied_filter_signal``."""
        if len(name) > 30:
            name = name[:30] + "..."
        self.filter_applied_friendly_name.setText(name)
        self.filter_applied_friendly_name.setEnabled(True)

    def _show_default_applied(self) -> None:
        """Show that no filter is applied. Use ``show_default_applied_signal``."""
        self.filter_applied_friendly_name.setText(self.DEFAULT_APPLIED_TEXT)
        self.filter_applied_friendly_name.setEnabled(False)

    def _show_validation_errors(self, errors: list[str]) -> None:
        """Show why the filter could not be applied, or hide the message if ``errors`` is empty."""
        self.filter_validation_errors.setText("\n".join(errors))
        self.filter_validation_errors.setVisible(bool(errors))
