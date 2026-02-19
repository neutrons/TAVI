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

   from tavi.meta.event.event_interface import Event
   from myapp.events import UserCreatedEvent
   from myapp.event_broker import EventBroker

   def on_user_created(event: UserCreatedEvent) -> None:
       print(f"User created: {event.user_id}")

   broker = EventBroker()
   broker.register(UserCreatedEvent, on_user_created)

Publishing Events
~~~~~~~~~~~~~~~~~

When an event is published, all subscribers registered for that event type are invoked.

.. code-block:: python

   event = UserCreatedEvent(user_id="123")
   broker.publish(event)

Each subscriber receives a **deep copy** of the event instance. This prevents
subscribers from mutating shared state and affecting other listeners.

Event Dispatch Semantics
------------------------

- **Dispatch is synchronous**: subscribers are called in the order they were registered.
- **Dispatch is type-based**: only subscribers registered for the exact event class
  (``type(event)``) are invoked.
- **Event instances are copied**: each subscriber receives an isolated event object.

Recursion Guard
---------------

The broker enforces a maximum call depth to prevent infinite or runaway recursion
when events trigger other events during handling.

If the maximum depth is exceeded, a ``RuntimeError`` is raised:

.. code-block:: python

   RuntimeError: Event recursive depth of 1 has been exceeded.

This protects against patterns like:

- A handler publishing the same event type it is subscribed to
- Circular event chains between handlers

If deeper event chaining is required, the maximum depth can be increased:

.. code-block:: python

   broker = EventBroker()
   broker.call_depth_max = 3

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
