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
       results = self.target(*self.args, **self.kwargs)
   except RecoverableError as e:
       ...
   except NonRecoverableError as e:
       ...
   except Exception as e:
       stack_trace = ...
       self.event_broker.publish(
           ExceptionEvent(e=NonRecoverableError(error_message, stack_trace))
       )

Key responsibilities of Worker:

- Prevent exceptions from escaping the thread
- Capture stack trace context
- Convert unexpected exceptions into domain-level ``TaviError``
- Dispatch via ``EventBroker``

Important:

If backend code already raises a ``TaviError`` subtype,
the Worker may pass it through directly instead of wrapping it.

Step 4 — Register a Handler
---------------------------

Handlers currently are registered in the frontend orchestration layer:

.. code-block:: python

   recovery.register(InvalidConfigurationError, self.handle_invalid_config)

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
       self.application_model.write_error_log(ex.stack_trace)
       TaviMessageBox.critical(self.view, "Error", str(ex))

This stops the workflow and returns control to the user.

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
- All errors become ``TaviError`` instances
- All errors enter the system through ``ExceptionEvent``
- Backend never interacts directly with UI

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
