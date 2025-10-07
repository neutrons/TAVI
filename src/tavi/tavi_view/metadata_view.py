from typing import Optional
from typing import Optional

from qtpy.QtCore import QObject
from qtpy.QtCore import QObject
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)


class MetaDataView(QWidget):
    """Main widget"""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the main widget
        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        self.display_metadata_callback = None

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.metadata_widget = MetaDataWidget(self)
        layout.addWidget(self.metadata_widget)

        # self.tree_widget = TreeViewWidget(self)
        # layout.addWidget(self.tree_widget)

        # self.load_widget.data_dir_or_files_signal.connect(self.load_view_data)



class MetaDataWidget(QWidget):
    """Widget that displays the metadata"""


    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the plotting widget

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutTop = QHBoxLayout()
        self.setLayout(layoutTop)

        self.filename_label = QLabel("File name:", self)
        self.filename_edit = QLineEdit(self)
        self.filename_edit.setEnabled(False)
        self.filename_label.setBuddy(self.filename_edit)

        layoutTop.addWidget(self.filename_label)
        layoutTop.addWidget(self.filename_edit)

    def set_values(self, values: str) -> None:
        self.filename_edit.setText(values)
