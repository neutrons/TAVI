EventBroker
===========

The ``EventBroker`` provides a lightweight publish/subscribe (pub-sub) mechanism for
decoupled event-driven communication between components. It acts as a central
dispatcher that routes events to registered subscribers based on event type.

Overview
--------

The broker is implemented as a singleton, so all parts of the application interact
with the same event registry. Components can:

- **Register** handlers for specific event types
- **Publish** events to notify all interested subscribers
- Rely on a simple recursion guard to prevent runaway event loops

This pattern is useful for cross-cutting concerns such as logging, state updates,
UI notifications, or domain events.

Basic Usage
-----------

Registering Subscribers
~~~~~~~~~~~~~~~~~~~~~~~

Subscribers are callables that accept a single event instance.

.. code-block:: python

   from tavi.meta.event.event_broker import EventBroker
   from tavi.meta.event.event_interface import Event

   class UserCreatedEvent(Event):
       user_id: str

   def on_user_created(event: UserCreatedEvent) -> None:
       print(f"User created: {event.user_id}")

   broker = EventBroker()
   broker.register(UserCreatedEvent, on_user_created)

Every event type must subclass :class:`tavi.meta.event.event_interface.Event`,
which is a pydantic ``BaseModel`` configured with
``arbitrary_types_allowed=True``. That is what makes the deep copy below
possible. TAVI's own event types live in ``tavi.meta.event.type``.

Publishing Events
~~~~~~~~~~~~~~~~~

When an event is published, all subscribers registered for that event type are invoked.

.. code-block:: python

   event = UserCreatedEvent(user_id="123")
   broker.publish(event)

Each subscriber receives a **deep copy** of the event instance
(``event.model_copy(deep=True)``). This prevents subscribers from mutating shared
state and affecting other listeners, and it is the mechanism that lets an event
safely carry model-owned objects — see
:doc:`../design/frontend/plot_data_model`.

Event Dispatch Semantics
------------------------

- **Dispatch is synchronous**: subscribers are called in the order they were registered.
- **Dispatch is type-based**: only subscribers registered for the exact event class
  (``type(event)``) are invoked.
- **Event instances are copied**: each subscriber receives an isolated event object.

Recursion Guard
---------------

The broker enforces a maximum call depth to prevent infinite or runaway recursion
when events trigger other events during handling. The default
``call_depth_max`` is **5**.

If the maximum depth is exceeded, ``publish`` raises:

.. code-block:: text

   RuntimeError: Event recursive depth of 5 has been exceeded.

This protects against patterns like:

- A handler publishing the same event type it is subscribed to
- Circular event chains between handlers

Note that legitimate chains count against this budget. The selection →
visualization flow already nests four deep on a fresh selection
(``FocusEvent`` → ``RawScanFocusEvent`` → ``PlotFocusEvent`` →
``ActivePlotChangedEvent``; see
:doc:`../design/frontend/visualization_flow`) — the limit was raised from 3
to 5 specifically to give that chain headroom once ``ActivePlotChangedEvent``
was added on the end of it. Switching the active plot from the "Current
Plot" dropdown is a separate, shallower chain (``FocusActivePlotEvent`` →
``ActivePlotChangedEvent``, 2 deep) that shares the same budget.

If deeper event chaining is required, raise the maximum on the singleton before
the chain runs:

.. code-block:: python

   broker = EventBroker()
   broker.call_depth_max = 5

Recommended Practices
---------------------

- **Keep handlers small and side-effect focused**
  Event handlers should perform limited, well-defined actions and avoid complex control flow.

- **Avoid cyclic event dependencies**
  Design event flows to be acyclic where possible. The recursion guard is a safety net,
  not a control mechanism.

- **Prefer domain-specific events**
  Use narrowly scoped event types (e.g., ``UserCreatedEvent`` instead of a generic
  ``UserEvent``) to keep subscriptions explicit and predictable.

- **Do not mutate incoming events**
  Although handlers receive copies, treat events as immutable to preserve intent
  and make behavior easier to reason about.

Typical Use Cases
-----------------

- Emitting domain events from application services
- Triggering side effects such as logging, metrics, or notifications
- Decoupling UI updates from core business logic
- Broadcasting lifecycle events (startup, shutdown, state changes)

Limitations
-----------

- No built-in support for asynchronous handlers
- No wildcard or base-class subscriptions (exact type matching only)
- No unregistration mechanism for subscribers
- No event classification system.  It does not validate that a subscriber *should* receive a specific event class. (Model vs Presenter)
- Global singleton scope may be undesirable in some testing or multi-tenant contexts

For more complex workflows (async dispatch, filtering, prioritization, or scoped
brokers), consider layering a more advanced event bus on top of this interface.
