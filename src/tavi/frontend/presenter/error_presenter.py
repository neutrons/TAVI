"""Presenter the orchestrates what frontend components to show when an error occurs."""

from tavi.backend.model.interface.application_model_interface import ApplicationModelInterface
from tavi.frontend.view.error_view import ErrorView
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError
from tavi.meta.exception.recovery_service import RecoveryService


class ErrorPresenter:
    """Presenter the orchestrates what frontend components to show when an error occurs."""

    def __init__(self, application_model: ApplicationModelInterface) -> None:
        """Initialise the error presenter."""
        super().__init__()
        self.recovery_service: RecoveryService = RecoveryService()
        self.application_model: ApplicationModelInterface = application_model

        self.recovery_service.register(NonRecoverableError, self.handle_nonrecoverable_exception)
        self.view = ErrorView()

    def handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
        """Emit signal to handle non-recoverable exception."""
        self.application_model.write_error_log(f"{ex.stack_trace}\n{str(ex)}")
        self.view.handle_nonrecoverable_exception(ex)
