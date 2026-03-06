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
- ``ExceptionEvent`` - Event wrapper for exceptions
- ``RecoveryService`` - Central exception router
- ``ErrorPresenter`` - Frontend error orchestration
- ``TaviMessageBox`` - UI dialog wrapper

Flow:

Exception → ``ExceptionEvent`` → ``RecoveryService`` → handler → frontend/backend action

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

- Subscribes to ``ExceptionEvent``
- Enforces ``TaviError``-only participation
- Routes exceptions by exact type
- Fails fast when no handler is registered

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

Currently, the presenter registers handling for ``NonRecoverableError`` and:

- Logs the error
- Displays a critical dialog
- Stops the current workflow


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
