Implementing a New Exception and Handler
========================================

This guide describes how to correctly introduce a new exception into the
backend → worker → recovery → handler pipeline.

The critical architectural rule is:

**Backend code raises exceptions.
Worker catches them.
RecoveryService routes them.
Handler defines policy.**

Exception Flow Overview
-----------------------

Backend Layer
    Raises domain exception.

Worker Thread Boundary
    Catches exception, wraps it as a ``TaviError`` (if needed),
    publishes ``ExceptionEvent``.

RecoveryService
    Routes exception to registered handler.

Frontend (ErrorPresenter)
    Defines user interaction and recovery policy.

Step 1 — Define the Exception
-----------------------------

Classify by workflow semantics:

- ``RecoverableError`` → additional system workflow required.
- ``NonRecoverableError`` → current workflow stops.

Example:

.. code-block:: python

   from tavi.meta.exception.nonrecoverable.base import NonRecoverableError

   class InvalidConfigurationError(NonRecoverableError):
       pass

Or:

.. code-block:: python

   from tavi.meta.exception.recoverable.base import RecoverableError

   class AuthenticationExpiredError(RecoverableError):
       pass

Exception classes should contain no behavior — only state.

Step 2 — Raise the Exception in Backend Code
--------------------------------------------

Backend services raise domain exceptions normally:

.. code-block:: python

   if not config.is_valid():
       raise InvalidConfigurationError(
           message="Configuration is invalid.",
           stack_trace="captured upstream or placeholder"
       )

No UI logic should exist in backend code.

Step 3 — Worker Catches and Dispatches
---------------------------------------

The ``Worker`` acts as the thread boundary and central exception capture point.

Simplified flow:

.. code-block:: python

   try:
       # Expects the return to be wrapped in a ModelResponse
       results: ModelResponse = self.target(*self.args, **self.kwargs)
   except Exception as e:  # noqa: BLE001
       stack_trace = f"{self.call_stack} \n {traceback.format_exc()}"
       error_message = str(e)
       results = ModelResponse(code=ResponseCode.ERROR, message=error_message)
       self.event_broker.publish(
           ExceptionEvent(error=NonRecoverableError(error_message, stack_trace))
       )

Key responsibilities of Worker:

- Prevent exceptions from escaping the thread
- Capture stack trace context — including the *submitting* thread's stack, saved
  as ``call_stack`` when ``WorkerPool.create_worker`` was called, which a bare
  ``format_exc()`` on the worker thread would not show
- Convert exceptions into domain-level ``TaviError``
- Dispatch via ``EventBroker``

.. important::

   The ``Worker`` currently has a **single** ``except Exception`` clause and
   wraps everything — including exceptions that are already ``TaviError``
   subtypes — in a fresh ``NonRecoverableError`` carrying only the original's
   ``str()``. There is no ``except RecoverableError`` branch and no pass-through
   for existing ``TaviError`` instances.

   Consequences to plan around:

   - A ``RecoverableError`` raised in backend code reaches the
     ``RecoveryService`` as a ``NonRecoverableError``, so a handler registered
     for the recoverable type will not fire.
   - A custom ``NonRecoverableError`` subtype loses its identity, so exact-type
     handler lookup falls back to the ``NonRecoverableError`` handler.

   Custom exception types are therefore only distinguishable today if they are
   raised and published on the caller's own thread rather than through a worker.
   Wiring per-type dispatch through the worker boundary requires teaching
   ``Worker.run`` to re-publish ``TaviError`` instances unchanged.

Note also that ``ExceptionEvent``'s field is named ``error``, not ``e``.

Step 4 — Register a Handler
---------------------------

Handlers currently are registered in the frontend orchestration layer
(``ErrorPresenter.__init__``):

.. code-block:: python

   self.recovery_service = RecoveryService()
   self.recovery_service.register(InvalidConfigurationError, self.handle_invalid_config)

``register`` raises ``RuntimeError`` if the type is not a ``TaviError``
subclass.

Rules:

- Exact type matching
- One handler per exception type
- Handler defines outcome (fatal, retry, ignore, escalate)

Step 5 — Implement the Handler
------------------------------

Handlers define policy and UI behavior.

Example: NonRecoverable

.. code-block:: python

   def handle_invalid_config(self, ex: InvalidConfigurationError) -> None:
       self.application_model.write_error_log(f"{ex.stack_trace}\n{str(ex)}")
       self._view.handle_invalid_config(ex)

This stops the workflow and returns control to the user.

.. warning::

   Do **not** call ``TaviMessageBox`` (or any other Qt widget) directly from a
   handler. Handlers run on whichever thread published the ``ExceptionEvent`` —
   normally a worker thread — and Qt widgets may only be touched from the GUI
   thread.

   ``ErrorPresenter`` delegates to ``ErrorView`` instead, which re-emits the
   exception as a Qt signal. The connected slot runs on the GUI thread and is
   the only place the message box is actually constructed:

   .. code-block:: python

      class ErrorView(QWidget):
          nonrecoverable_signal = Signal(NonRecoverableError)

          def __init__(self) -> None:
              super().__init__()
              self.nonrecoverable_signal.connect(self._handle_nonrecoverable_exception)

          def handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
              self.nonrecoverable_signal.emit(ex)          # any thread

          def _handle_nonrecoverable_exception(self, ex: NonRecoverableError) -> None:
              TaviMessageBox.critical(self, "Error", str(ex))   # GUI thread

Example: Recoverable

.. code-block:: python

   def handle_auth_expired(self, ex: AuthenticationExpiredError) -> None:
       if self.auth_service.refresh_token():
           self.retry_original_action()
       else:
           self.escalate_to_fatal(ex)

Worker Contract
---------------

The Worker enforces:

- All backend exceptions are captured
- All errors become ``NonRecoverableError`` instances (see the caveat in Step 3)
- All errors enter the system through ``ExceptionEvent``
- Backend never interacts directly with UI

After publishing, the worker still emits ``finished`` and then asserts that the
target returned a ``ModelResponse``. On the error path it substitutes
``ModelResponse(code=ResponseCode.ERROR, message=...)`` itself, so the assertion
only fires for a target that returned successfully with the wrong type.

This keeps thread boundaries clean and recovery centralized.

Design Guidelines
-----------------

1. Raise exceptions only in backend logic.
2. Never perform UI operations inside backend services.
3. Let Worker normalize and dispatch errors.
4. Keep exception classes minimal.
5. Keep handlers focused and deterministic.
6. Always register new exception types explicitly.

Checklist
---------

When adding a new exception:

- [ ] Subclass correct base (Recoverable vs NonRecoverable)
- [ ] Raise only in backend layer
- [ ] Ensure Worker dispatches it correctly
- [ ] Register handler in presenter
- [ ] Define logging policy
- [ ] Define escalation behavior

Architectural Summary
---------------------

Backend expresses failure.
Worker captures and dispatches.
RecoveryService routes.
Handler decides outcome.

This keeps responsibilities sharply separated and prevents cross-layer leakage.
