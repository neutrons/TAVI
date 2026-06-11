"""Data file view widget."""

from typing import Any

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QAbstractButton,
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

    def __init__(self, parent: Any = None) -> None:
        """Construct data file view."""
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the data file UI."""
        layout = QVBoxLayout(self)

        top_split = QSplitter(Qt.Orientation.Horizontal)

        # Data table
        table = QTableWidget(10, 5)
        table.setAlternatingRowColors(True)
        table.setHorizontalHeaderLabels(["detector", "s1", "s2", "def_x", "def_y"])
        table.setCornerButtonEnabled(True)

        top_split.addWidget(table)

        # Variables list
        variable_table = QTableWidget()
        variable_table.setColumnCount(1)
        variable_table.setAlternatingRowColors(True)
        variable_table.setHorizontalHeaderLabels(["Variables"])
        variable_table.setCornerButtonEnabled(True)
        variable_table.verticalHeader().setMinimumWidth(20)
        for v in ["detector", "s1", "s2", "sample", "monitor", "mcu"]:
            item = self._add_row(variable_table, v)
            item.setCheckState(Qt.CheckState.Checked)
            item.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            item.setCheckState(Qt.CheckState.Unchecked)

        button = variable_table.findChild(QAbstractButton)
        if button:
            checkbox = QCheckBox("", button)
            checkbox.setMinimumSize(button.width() / 2, button.height())

            def flipChecks(checkState: Any) -> None:
                for row in range(variable_table.rowCount()):
                    variable_table.item(row, 0).setCheckState(checkState)

            checkbox.checkStateChanged.connect(flipChecks)

        top_split.addWidget(variable_table)
        layout.addWidget(top_split, 2)

        # Metadata
        meta_tabs = QTabWidget()
        meta_tabs.addTab(self._build_scan_info_tab(), "Scan Info")
        layout.addWidget(meta_tabs, 1)

        # Buttons
        btn_layout = QHBoxLayout()
        btn_layout.addWidget(QPushButton("Restore Metadata From File"))
        btn_layout.addWidget(QPushButton("Save Modified Metadata To File"))
        layout.addLayout(btn_layout)

    def _add_row(self, table: QTableWidget, *args: Any) -> QTableWidgetItem:
        """Add a row to a table."""
        row_idx = table.rowCount()
        table.insertRow(row_idx)
        item = None
        for col_idx, value in enumerate(args):
            item = QTableWidgetItem(str(value))
            table.setItem(row_idx, col_idx, item)
        return item

    def _build_scan_info_tab(self) -> QTableWidget:
        """Build the scan info metadata tab."""
        metadata = QTableWidget()
        metadata.setColumnCount(2)
        metadata.setAlternatingRowColors(True)

        self._add_row(metadata, "scan number", "1")
        self._add_row(metadata, "scan title", "title")
        self._add_row(metadata, "preset channel", "MCU")
        self._add_row(metadata, "preset type", "normal")
        self._add_row(metadata, "preset value", "10")
        self._add_row(metadata, "def_x", "s2")
        self._add_row(metadata, "def_y", "detector")
        self._add_row(metadata, "start time", "")
        self._add_row(metadata, "end time", "")
        self._add_row(metadata, "date", "")

        for row in range(metadata.rowCount()):
            item = metadata.item(row, 0)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)

        header = metadata.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        header.hide()

        metadata.verticalHeader().hide()

        return metadata
