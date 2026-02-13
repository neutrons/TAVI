"""Presenter the orchestrates what frontend components to show when an error occurs."""

from qtpy.QtCore import QObject, Signal

from tavi.frontend.view.error_view import ErrorView
from tavi.frontend.widget.tavi_message_box import TaviMessageBox
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError
from tavi.meta.exception.recovery_service import RecoveryService


class ErrorPresenter(QObject):
    """Presenter the orchestrates what frontend components to show when an error occurs."""

    nonrecoverable_signal = Signal(NonRecoverableError)

    def __init__(self) -> None:
        """Initialise the error presenter."""
        super().__init__()
        self.recovery_service = RecoveryService()

        self.nonrecoverable_signal.connect(self._handle_nonrecoverable_exception)

        self.recovery_service.register(NonRecoverableError, self.handle_nonrecoverable_exception)
        self.view = ErrorView()

    def handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
        """Emit signal to handle non-recoverable exception."""
        self.nonrecoverable_signal.emit(ex)

    def _handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
        """Handle non-recoverable exception."""
        TaviMessageBox.critical(self.view, "Error", str(ex))
