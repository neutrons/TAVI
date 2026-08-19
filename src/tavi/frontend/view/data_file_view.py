"""Data file view widget."""

from typing import Any

from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QGuiApplication, QKeySequence
from qtpy.QtWidgets import (
    QAbstractButton,
    QAction,
    QCheckBox,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)


class DataFileView(QWidget):
    """Data file panel widget."""

    title_changed = Signal(str)
    scan_focus_changed = Signal(object)

    def __init__(self, parent: Any = None) -> None:
        """Construct data file view."""
        super().__init__(parent)
        self._build_ui()
        # AutoConnection: direct call on the GUI thread (tests), queued hop when emitted from
        # a worker thread (PlotModel/TaviProjectModel running behind their Proxy).
        self.scan_focus_changed.connect(self._apply_scan_focus)

    def _build_ui(self) -> None:
        """Build the data file UI."""
        layout = QVBoxLayout(self)

        top_split = QSplitter(Qt.Orientation.Horizontal)

        # Data table
        self.data_table = QTableWidget(0, 0)
        self.data_table.setAlternatingRowColors(True)
        self.data_table.setCornerButtonEnabled(True)
        self.data_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        self._enable_copy(self.data_table)

        top_split.addWidget(self.data_table)

        # Variables list
        self.variable_table = QTableWidget()
        self.variable_table.setColumnCount(1)
        self.variable_table.setAlternatingRowColors(True)
        self.variable_table.setHorizontalHeaderLabels(["Variables"])
        self.variable_table.setCornerButtonEnabled(True)
        self.variable_table.verticalHeader().setMinimumWidth(20)
        self.variable_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)

        button = self.variable_table.findChild(QAbstractButton)
        if button:
            checkbox = QCheckBox("", button)
            checkbox.setMinimumSize(button.width() // 2, button.height())

            def flipChecks(checkState: Any) -> None:
                for row in range(self.variable_table.rowCount()):
                    self.variable_table.item(row, 0).setCheckState(checkState)

            checkbox.stateChanged.connect(flipChecks)

        self.variable_table.itemChanged.connect(self._on_variable_check_changed)
        self.variable_table.verticalHeader().setSectionsMovable(True)
        self.variable_table.verticalHeader().sectionMoved.connect(self._on_variable_row_moved)

        top_split.addWidget(self.variable_table)
        top_split.setStretchFactor(0, 2)
        top_split.setStretchFactor(1, 1)
        layout.addWidget(top_split, 2)

        # Metadata
        self.meta_tabs = QTabWidget()
        self._show_empty_metadata_tab()
        layout.addWidget(self.meta_tabs, 1)

        # Buttons
        btn_layout = QHBoxLayout()
        self.restore_metadata_button = QPushButton("Restore Metadata From File")
        self.restore_metadata_button.setEnabled(False)
        self.save_metadata_button = QPushButton("Save Modified Metadata To File")
        self.save_metadata_button.setEnabled(False)
        btn_layout.addWidget(self.restore_metadata_button)
        btn_layout.addWidget(self.save_metadata_button)
        layout.addLayout(btn_layout)

    def _enable_copy(self, table: QTableWidget) -> None:
        """Give a table a Ctrl+C / right-click "Copy" action - item views have no copy support of their own."""
        copy_action = QAction("Copy", table)
        copy_action.setShortcut(QKeySequence.StandardKey.Copy)
        copy_action.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
        copy_action.triggered.connect(lambda: self._copy_selection(table))
        table.addAction(copy_action)
        table.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)

    def _copy_selection(self, table: QTableWidget) -> None:
        """
        Put the selected cells on the clipboard as tab-separated text, one line per row.

        Rows and columns are emitted in the header's *visual* order so a copy matches what
        the user sees after dragging variables around, and hidden columns are left out.
        Cells inside the selection's bounding box that are not themselves selected come
        out empty, which keeps non-rectangular selections aligned.
        """
        selection = table.selectedRanges()
        if not selection:
            return

        rows = sorted({row for rng in selection for row in range(rng.topRow(), rng.bottomRow() + 1)})
        columns = {col for rng in selection for col in range(rng.leftColumn(), rng.rightColumn() + 1)}
        columns = sorted(
            (col for col in columns if not table.isColumnHidden(col)),
            key=table.horizontalHeader().visualIndex,
        )

        lines = []
        for row in rows:
            if table.isRowHidden(row):
                continue
            cells = []
            for col in columns:
                item = table.item(row, col)
                cells.append(item.text() if item is not None and item.isSelected() else "")
            lines.append("\t".join(cells))

        QGuiApplication.clipboard().setText("\n".join(lines))

    def set_title(self, title: str) -> None:
        """Emit ``title_changed`` so an owning tab widget can relabel itself."""
        self.title_changed.emit(title)

    def _apply_scan_focus(self, scan: Any) -> None:
        """Repopulate every data widget section from a newly-focused scan, or clear them if there isn't one."""
        if scan is None:
            self.clear_data()
            self.set_title("Data File")
            return
        self.populate_columns(scan.data.data)
        self.populate_variables(list(scan.data.data.keys()))
        self.populate_metadata(scan.metadata.by_category())
        self.set_title(f"Data File ({scan.tavimeta.friendly_name})")

    def populate_columns(self, data: dict[str, list[float]]) -> None:
        """Repopulate the data table with a newly-focused scan's column values."""
        columns = list(data.keys())
        self.data_table.setRowCount(0)
        self.data_table.setColumnCount(len(columns))
        self.data_table.setHorizontalHeaderLabels(columns)
        row_count = max((len(values) for values in data.values()), default=0)
        self.data_table.setRowCount(row_count)
        for col_idx, name in enumerate(columns):
            for row_idx, value in enumerate(data[name]):
                item = QTableWidgetItem(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.data_table.setItem(row_idx, col_idx, item)

    def populate_variables(self, names: list[str]) -> None:
        """Repopulate the variable list with a newly-focused scan's column names."""
        self.variable_table.setRowCount(0)
        for name in names:
            item = self._add_row(self.variable_table, name)
            item.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            item.setCheckState(Qt.CheckState.Checked)

    def _on_variable_check_changed(self, item: QTableWidgetItem) -> None:
        """Show/hide the data table column matching a toggled variable's name."""
        visible = item.checkState() == Qt.CheckState.Checked
        for col in range(self.data_table.columnCount()):
            header_item = self.data_table.horizontalHeaderItem(col)
            if header_item is not None and header_item.text() == item.text():
                self.data_table.setColumnHidden(col, not visible)
                break

    def _on_variable_row_moved(self, logical_index: int, old_visual_index: int, new_visual_index: int) -> None:
        """
        Mirror a dragged variable row onto the matching data table column.

        ``variable_table`` rows and ``data_table`` columns are always populated together
        from the same ordered name list, so a row's logical index is always the same
        column's logical index — the stable identity used to find it regardless of
        either header's current visual order.
        """
        data_header = self.data_table.horizontalHeader()
        current_visual = data_header.visualIndex(logical_index)
        data_header.moveSection(current_visual, new_visual_index)

    def populate_metadata(self, metadata: dict[str, Any]) -> None:
        """Rebuild the metadata tab widget: one tab per top-level key, one table row per k/v pair."""
        self.meta_tabs.clear()
        for tab_name, fields in metadata.items():
            if not isinstance(fields, dict):
                raise ValueError(f"Metadata tab {tab_name!r} must be a dict, not {type(fields).__name__!r}.")
            self.meta_tabs.addTab(self._build_metadata_tab(fields), tab_name)

    def clear_data(self) -> None:
        """Clear the data table, variable list, and metadata tabs, e.g. when a focus event carries no scans."""
        self.data_table.setRowCount(0)
        self.data_table.setColumnCount(0)
        self.variable_table.setRowCount(0)
        self._show_empty_metadata_tab()

    def _show_empty_metadata_tab(self) -> None:
        """Reset the metadata tab widget to a single empty 'Empty' tab."""
        self.meta_tabs.clear()
        empty_table = QTableWidget()
        self._finalize_metadata_table(empty_table)
        self.meta_tabs.addTab(empty_table, "Empty")

    def _add_row(self, table: QTableWidget, *args: Any) -> QTableWidgetItem:
        """Add a row to a table. Entries are non-editable for now."""
        row_idx = table.rowCount()
        table.insertRow(row_idx)
        item = None
        for col_idx, value in enumerate(args):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_idx, col_idx, item)
        return item

    def _build_metadata_tab(self, fields: Any) -> QTableWidget:
        """Build one metadata tab's table - a k/v mapping renders 2-column, a bare list renders 1-column."""
        if isinstance(fields, dict):
            return self._build_metadata_mapping_table(fields)
        return self._build_metadata_list_table(fields)

    def _build_metadata_mapping_table(self, fields: dict[str, Any]) -> QTableWidget:
        """Build a 2-column key/value table."""
        table = QTableWidget()
        table.setColumnCount(2)
        table.setAlternatingRowColors(True)

        for key, value in fields.items():
            self._add_row(table, key, value)

        self._finalize_metadata_table(table)
        return table

    def _build_metadata_list_table(self, items: list) -> QTableWidget:
        """Build a 1-column table, one row per list entry."""
        table = QTableWidget()
        table.setColumnCount(1)
        table.setAlternatingRowColors(True)

        for value in items:
            self._add_row(table, value)

        self._finalize_metadata_table(table)
        return table

    def _finalize_metadata_table(self, table: QTableWidget) -> None:
        """Apply shared header/sizing styling to a metadata tab's table."""
        header = table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        header.hide()

        table.verticalHeader().hide()
        self._enable_copy(table)
