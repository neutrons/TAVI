"""A column-name field, shared by the panels that ask the user to name a scan column."""

from qtpy.QtWidgets import QComboBox


def typed_column_combo(default: str = "") -> QComboBox:
    """
    Build a column-name field: free text, with the focused scan's own columns offered as a dropdown.

    Typed rather than picked-only because a name can be set before the scan that has it is
    focused, and a scan need not advertise every column that can be asked for. NoInsert keeps
    what was typed out of the list, which describes the scan rather than what the user tried.
    """
    combo = QComboBox()
    combo.setEditable(True)
    combo.setInsertPolicy(QComboBox.NoInsert)
    combo.setEditText(default)
    return combo


def set_combo_text(combo: QComboBox, text: str) -> None:
    """
    Show ``text`` in an editable combo, keeping its index in step, and emitting nothing.

    ``setCurrentText`` only writes the line edit on an editable combo and leaves the index
    where it was, so the two drift apart - and picking the item that happens to sit at the
    stale index then changes no index, firing no replot.
    """
    combo.blockSignals(True)
    try:
        # -1 when the name is not one of the focused scan's columns, which is allowed here.
        combo.setCurrentIndex(combo.findText(text))
        combo.setEditText(text)
    finally:
        combo.blockSignals(False)


def set_combo_options(combo: QComboBox, columns: list[str]) -> None:
    """Replace a column field's dropdown with ``columns``, keeping whatever name it already shows."""
    current = combo.currentText()
    combo.blockSignals(True)
    try:
        combo.clear()
        combo.addItems(columns)
    finally:
        combo.blockSignals(False)
    set_combo_text(combo, current)
