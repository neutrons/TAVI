from typing import Optional

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)


class _UiBridge(QObject):
    """
    Thread-safe bridge to deliver updates on the GUI thread. Qt forbitds modifying
    UI elements from a different thread. The data needs to be passed as a signal.

    """

    set_metadata_signal = Signal(str)


class MetaDataView(QWidget):
    """
    Prototypical widget responsible for displaying metadata associated
    with loaded scans or files.

    Parameters
    ----------
    parent : Optional[QObject], default=None
        Parent widget used for Qt ownership and memory management.

    Attributes
    ----------
    metadata_widget : MetaDataWidget
        The inner widget responsible for displaying filename and related metadata.
    _bridge : _UiBridge
        Thread-safe bridge for routing metadata updates to the GUI.
    display_metadata_callback : Callable or None
        Placeholder for a presenter-provided callback (future use).
    """

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

        self._bridge = _UiBridge()
        self._bridge.set_metadata_signal.connect(
            self.metadata_widget.set_values,
            type=Qt.QueuedConnection,  # run safely on GUI thread
        )

    def update_metadata(self, event: None) -> None:
        """
        Request a metadata update by emitting the change through the `_UiBridge`.

        This method should be called by the presenter when the model produces
        new metadata.

        Parameters
        ----------
        event : str
            The metadata value to display. Typically the filename or metadata
            summary associated with the currently selected scan.
        """
        self._bridge.set_metadata_signal.emit(event)


class MetaDataWidget(QWidget):
    """
    Widget responsible for displaying metadata fields (e.g., filename).

    The widget contains a label and a disabled QLineEdit. The presenter or
    parent view updates its contents through `set_values()`, typically routed
    via `_UiBridge` to ensure thread safety.
    """

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the metadata widget

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutTop = QHBoxLayout()
        self.setLayout(layoutTop)

        self.filename_label = QLabel("File name:", self)
        self.filename_edit = QLineEdit(self)
        self.filename_edit.setStyleSheet("color: black;")
        self.filename_edit.setEnabled(False)
        self.filename_label.setBuddy(self.filename_edit)

        layoutTop.addWidget(self.filename_label)
        layoutTop.addWidget(self.filename_edit)

    def set_values(self, values: str) -> None:
        """
        Update the metadata display field.

        Parameters
        ----------
        values : str
            Text to display in the filename field (e.g., the selected scan name).
        """
        self.filename_edit.setText(values)
