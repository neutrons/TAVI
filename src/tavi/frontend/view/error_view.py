"""View to show when an error occurs."""

from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import QWidget

from tavi.frontend.widget.tavi_message_box import TaviMessageBox
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError


class ErrorView(QWidget):
    """Placeholder view for when an error occurs."""

    nonrecoverable_signal = Signal(NonRecoverableError)

    def __init__(self) -> None:
        """Initialize error handling and remove clickabilitiy."""
        super().__init__()
        # PyQt / PySide
        self.setAttribute(Qt.WA_TransparentForMouseEvents)

        self.nonrecoverable_signal.connect(self._handle_nonrecoverable_exception)

    def handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
        """Emit signal to handle non-recoverable exception."""
        self.nonrecoverable_signal.emit(ex)

    def _handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
        """Handle non-recoverable exception."""
        TaviMessageBox.critical(self, "Error", str(ex))
