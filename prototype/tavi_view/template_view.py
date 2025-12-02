from typing import Optional

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from tavi.EventBroker.event_type import template_data


class _UiBridge(QObject):
    """
    Thread-safe bridge to deliver updates on the GUI thread. Qt forbitds modifying
    UI elements from a different thread. The data needs to be passed as a signal.

    """

    set_template_signal = Signal(str)


class TemplateView(QWidget):
    """
    A template of a widget that will react to a secondary model of the TaviProjectModel.

    This view contains a `TemplateWidget` (which provides the actual display UI)
    and a `_UiBridge` to ensure updates from presenters or background operations
    are applied safely on the GUI thread.

    Attributes
    ----------
    template_widget : TemplateWidget
        The inner display widget showing the template value.
    _bridge : _UiBridge
        Thread-safe event bridge used to route updates to the GUI.
    """

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """
        Initialize the template view, set up layouts, and connect the bridge
        signal to the template widget’s update method.
        """
        super().__init__(parent)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.template_widget = TemplateWidget(self)
        layout.addWidget(self.template_widget)

        self._bridge = _UiBridge()
        self._bridge.set_template_signal.connect(
            self.template_widget.set_values,
            type=Qt.QueuedConnection,
        )

    def update_template(self, event: template_data) -> None:
        """
        Request an update to the displayed template value.

        This method should be called by a presenter when new template information
        becomes available. It emits the update through `_UiBridge`, ensuring the
        update executes on the GUI thread.

        Parameters
        ----------
        event : str
            The new template text to display.
        """
        self._bridge.set_template_signal.emit(event)


class TemplateWidget(QWidget):
    """
    Widget responsible for displaying a single template field (e.g., the
    “next file name” that the system expects or will generate).

    The widget consists of a label and a read-only text field that displays the
    current template value.
    """

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the template widget

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutTop = QHBoxLayout()
        self.setLayout(layoutTop)

        self.filename_label = QLabel("Next file:", self)
        self.filename_edit = QLineEdit(self)
        self.filename_edit.setStyleSheet("color: black;")
        self.filename_edit.setEnabled(False)
        self.filename_label.setBuddy(self.filename_edit)

        layoutTop.addWidget(self.filename_label)
        layoutTop.addWidget(self.filename_edit)

    def set_values(self, values: str) -> None:
        """
        Update the displayed template value.

        Parameters
        ----------
        values : str
            The template text to display (typically a filename or prefix).
        """
        self.filename_edit.setText(values)
