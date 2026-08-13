Error Recovery System
=====================

The error recovery system provides a centralized, event-driven mechanism for
handling domain exceptions while cleanly separating backend logic from frontend
orchestration.

Exceptions are not propagated directly to the UI. Instead, they are:

1. Modeled as domain exceptions
2. Published as events
3. Routed through a recovery service
4. Interpreted and acted upon by discrete handlers

This architecture allows recovery *policy* to remain flexible and contextual.

Core Concepts
-------------

This system distinguishes **workflow semantics** from **severity policy**.

Workflow Classification
~~~~~~~~~~~~~~~~~~~~~~~

Exceptions are classified by whether the system can continue pursuing the
current user goal.

- ``RecoverableError``
  The application must execute additional workflow in order to continue
  toward the same user goal.

- ``NonRecoverableError``
  The current workflow stops immediately. The user must re-initiate the action,
  even if the error itself is user-correctable.

Architecture Overview
---------------------

Components:

- ``TaviError`` - Base exception for all domain errors
- ``RecoverableError`` - Requires additional system workflow
- ``NonRecoverableError`` - Stops the current workflow
- ``ExceptionEvent`` - Event wrapper for exceptions (field: ``error``)
- ``Worker`` - Thread boundary; captures every backend exception and publishes it
- ``RecoveryService`` - Central exception router
- ``ErrorPresenter`` - Frontend error orchestration
- ``ErrorView`` - Marshals the handler back onto the GUI thread via a Qt signal
- ``TaviMessageBox`` - UI dialog wrapper

Flow:

Exception → ``Worker`` → ``ExceptionEvent`` → ``RecoveryService`` → handler → ``ErrorPresenter`` → ``ErrorView`` → ``TaviMessageBox``

Exception Hierarchy
-------------------

TaviError
~~~~~~~~~

Base class for all framework-level exceptions.

Attributes:

- ``message`` - User-facing description
- ``stack_trace`` - Captured traceback

Implements ``__deepcopy__`` for safe event propagation.

All exceptions participating in recovery must inherit from ``TaviError``.

RecoverableError
~~~~~~~~~~~~~~~~

Represents an error where the application can continue toward the same user goal
by executing additional workflow.

Examples:

- Refreshing expired credentials
- Resolving a missing dependency
- Retrying a transient failure
- Requesting additional information required to proceed

NonRecoverableError
~~~~~~~~~~~~~~~~~~~

Represents an error that stops the current workflow immediately.

Examples:

- Validation failures
- Invalid user input
- Business rule violations
- Preconditions not met

These errors may still be fatal *depending on how they are handled*.

RecoveryService
---------------

The ``RecoveryService`` is responsible for routing exceptions to handlers.

Responsibilities:

- Subscribes to ``ExceptionEvent`` on the ``EventBroker``
- Enforces ``TaviError``-only participation at *registration* time —
  ``register`` raises ``RuntimeError`` for anything that is not a ``TaviError``
  subclass
- Routes exceptions by exact type (``type(ex)``, no inheritance traversal)
- Fails fast when no handler is registered: ``default_handler`` raises
  ``RuntimeError(f"FATAL: no handler for exception found {ex}")``

Registration
~~~~~~~~~~~~

.. code-block:: python

   recovery.register(MyErrorType, handler)

Constraints:

- One handler per exception type
- Exact type matching (no inheritance traversal)
- Handler decides severity and outcome

Dispatch Semantics
~~~~~~~~~~~~~~~~~~

When an exception event is published:

1. The exception is extracted
2. A handler is resolved
3. The handler executes
4. Outcome is entirely handler-defined

No assumptions are made by the recovery system about fatality, retries, or UI.

ErrorPresenter (Frontend Orchestration)
---------------------------------------

The ``ErrorPresenter`` is the frontend implementation of error handling policy.

Responsibilities:

- Registering handlers with ``RecoveryService``
- Translating exceptions into UI actions
- Logging errors
- Determining whether the application continues or halts

Current Behavior
~~~~~~~~~~~~~~~~

Currently, the presenter registers exactly one handler — for
``NonRecoverableError`` — which:

- Writes ``stack_trace`` plus the message to a timestamped error log via
  ``ApplicationModel.write_error_log``
- Forwards the exception to ``ErrorView.handle_nonrecoverable_exception``
- Stops the current workflow

``ErrorView`` does not show the dialog directly either. Because the exception
arrives on a worker thread, the view re-emits it as a Qt signal
(``nonrecoverable_signal``), and the connected slot — now running on the GUI
thread — calls ``TaviMessageBox.critical(self, "Error", str(ex))``.

No handler is registered for ``RecoverableError`` yet, so one reaching the
service today would hit ``default_handler`` and raise.


Future Recoverable Workflow Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A recoverable flow may:

1. Intercept a ``RecoverableError``
2. Execute additional orchestration
3. Resume the original operation
4. Escalate to fatal if recovery fails


Design Characteristics
----------------------

- Event-driven and synchronous
- Centralized policy definition
- Explicit handler registration
- Clear separation of domain errors and UX decisions

Recommended Practices
---------------------

- Classify errors by *workflow semantics*, not severity
- Treat fatality as a handler-level decision
- Keep handlers small and explicit
- Escalate intentionally, not implicitly
- Document recovery intent per exception type

This model preserves flexibility while keeping semantics precise:

**Types describe workflow impact.
Handlers decide severity and outcome.**
