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

   class AddUserEvent(Event):
       user_id: str

   def on_user_created(event: AddUserEvent) -> None:
       print(f"User added: {event.user_id}")

   broker = EventBroker()
   broker.register(AddUserEvent, on_user_created)

Every event type must subclass :class:`tavi.meta.event.event_interface.Event`,
which is a pydantic ``BaseModel`` configured with
``arbitrary_types_allowed=True``. That is what makes the deep copy below
possible. TAVI's own event types live in ``tavi.meta.event.type``.

Publishing Events
~~~~~~~~~~~~~~~~~

When an event is published, all subscribers registered for that event type are invoked.

.. code-block:: python

   event = AddUserEvent(user_id="123")
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

Note that legitimate chains count against this budget. The deepest chains in
TAVI are a fresh tree selection, which runs ``FocusEvent`` (1) →
``FocusRawScanEvent`` or ``FocusFitEvent`` (2) → ``FocusPlotEvent`` (3) →
``StageSeriesEvent`` (4) → ``SyncStageEvent`` (5); see
:doc:`../design/frontend/visualization_flow`. ``SyncStageEvent`` handlers
therefore run at the limit and must not publish anything themselves. Staging
from the "Current Plot" dropdown or the "Apply All" checkbox is a separate,
shallower chain (``ClearStageEvent`` and ``StageSeriesEvent`` (1) →
``SyncStageEvent`` (2)) that shares the same budget.

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
  Use narrowly scoped event types (e.g., ``AddUserEvent`` instead of a generic
  ``UserEvent``) to keep subscriptions explicit and predictable.

- **Follow the event lexicon**
  Name every event ``VerbNounEvent`` with a verb from `Event lexicon`_, keep each
  one a single prime operation, and use a direct model call rather than an event
  for a request only one model answers. `Designing TAVI events`_ explains why.

- **Do not mutate incoming events**
  Although handlers receive copies, treat events as immutable to preserve intent
  and make behavior easier to reason about.

Designing TAVI events
---------------------

TAVI's events are a small vocabulary. Every event names one operation on one
concept, so a chain of events reads like a sentence, and a subscriber can react
to exactly the operation it cares about without knowing who else is listening.

Naming
~~~~~~

Every event is named ``VerbNounEvent``: the operation first, then what it
applies to (``FocusPlotEvent``, ``SyncStageEvent``, ``RemoveFitEvent``). The
verb comes from the lexicon below. Names describe the domain, never a widget:
``ClearFocusEvent``, not ``ResetFitPanelEvent``, and ``StageSeriesEvent``, not
``ApplyAllChangedEvent``. A model publishing "reset the panel" would know a
panel exists; each presenter decides for itself what an operation means for
its own view.

Concepts
~~~~~~~~

Three nested sets describe what the user is working with:

**Focus**
    Everything the user is looking at: the scans, plots and fits selected in
    the project tree, and the series they put on the plotter. Focus is for
    looking.

**Stage**
    The part of the focus that edits and fits apply to. With "Apply All"
    checked every focused series is staged; without it, only the one picked in
    the "Current Plot" dropdown. Several staged series are fit in sequence.

**Lead**
    The first staged series. It is the one the plotter's fields, the data tab
    and the fitting panel show. There is no separate "active" series: with one
    series staged, the lead is that series.

``PlotModel`` owns the focused plots and the stage. The user changes both
through UI interaction, so presenters publish the requests, and ``PlotModel``
syncs the result back out.

Event lexicon
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 12 58 30

   * - Verb
     - Exact meaning
     - Events
   * - Focus
     - Add items to what the user is looking at. Focus events only ever add;
       they never replace what is already focused.
     - ``FocusEvent``, ``FocusRawScanEvent``, ``FocusPlotEvent``,
       ``FocusFitEvent``
   * - Stage
     - Add focused series to the stage. Like focus, staging only adds.
     - ``StageSeriesEvent``
   * - Clear
     - Empty a set so it can be rebuilt. Always paired with the additive verb
       for the same set. Clearing the focus also clears the stage, since the
       stage is part of the focus.
     - ``ClearFocusEvent``, ``ClearStageEvent``
   * - Sync
     - A model's answer: "this is how X stands now". Display only. Carries the
       data the UI needs, so presenters never reach into a model. Sent after
       the state changes, whatever changed it.
     - ``SyncPlotEvent``, ``SyncStageEvent``, ``SyncFitEvent``,
       ``SyncFitSpecEvent``, ``SyncFitHistoryEvent``, ``SyncPeakParamsEvent``,
       ``SyncBackgroundParamsEvent``, ``SyncRecentProjectsEvent``
   * - Recompute
     - Ask the model that computes X to compute it again for display. Never
       writes to the project.
     - ``RecomputeFitEvent``
   * - Save
     - Write X into the project (``TaviData``).
     - ``SavePlotEvent``, ``SaveFitEvent``
   * - Add / Remove
     - Announce that X has entered or left the project.
     - ``AddRawScanEvent``, ``AddPlotEvent``, ``AddFitEvent``,
       ``RemoveRawScanEvent``, ``RemovePlotEvent``, ``RemoveFitEvent``
   * - Restore
     - Announce that X was put back to an earlier saved state.
     - ``RestoreFitMemberEvent``
   * - Set
     - A display preference passed between presenters. No model state changes.
     - ``SetFitComponentsVisibleEvent``
   * - Report
     - Surface an error to the user.
     - ``ReportErrorEvent``
   * - Start
     - An application lifecycle step.
     - ``StartApplicationEvent``

``FocusEvent`` is the one generic event. The project tree publishes the bare
uuids the user selected, because only ``TaviData`` knows which type each uuid
is. ``TaviProjectModel`` narrows it into a ``ClearFocusEvent`` followed by one
specific focus event per type.

Prime events, not compound ones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each event is a prime: one operation that can't be split further. A user
action that is really several operations publishes several events, in order.
Changing the selection is a clear followed by a focus, so it is
``ClearFocusEvent`` and then ``FocusPlotEvent``, never one
``ChangeSelectionEvent``. Restaging is ``ClearStageEvent`` and then
``StageSeriesEvent``.

Compound events hide operations from subscribers that need to react to them
differently. Signs that an event is compound:

* A past-tense "changed" name (``ActivePlotChangedEvent``). "Changed" usually
  means "cleared, then set".
* A flag that changes what the event means (``exclusive`` meant "clear first").
* An optional payload where ``None`` means something different
  (``scan=None`` meant "nothing is active any more").
* A field that carries a second operation along (``also_plots`` focused plots
  inside a raw-scan focus, so the two wouldn't overwrite each other).
* One event for both a new state and a refresh of the same state. A new focus
  resets panels; a refresh after an edit must not. That is why a refresh is
  ``SyncPlotEvent`` and not another ``FocusPlotEvent``.

Making the additive verbs only add is what lets primes compose. Independent
publishers can each focus their own items after one clear, without
overwriting each other and without knowing about each other.

Keep a compound event only when the combination needs a reaction that no
sequence of its primes can produce. Wanting to render once rather than twice
is not such a reaction: it is an optimisation, not a different meaning.

Keep events few and high-level. Add a new event when a new concept or
operation appears, not to save a subscriber a line of code.

Events or model calls
~~~~~~~~~~~~~~~~~~~~~

Events are for broadcasting a concept to everyone who cares: one to many.
Focus, stage, sync, add and remove are concepts any number of presenters and
models may react to.

A request that exactly one model answers, because it changes state only that
model owns, is many to one. It is a direct call through that model's
``Proxy`` instead, such as ``PlotModel.update_fields``,
``PlotModel.set_show_title``, ``FitModel.perform_fit``,
``FitModel.sync_fit_spec`` and ``TaviProjectModel.undo_fit_member``. The
model then announces the result with a sync event, which every interested
subscriber receives. Broadcasting such a request would only add a hop, and
would hide who is responsible for answering it.

Staging is the exception that proves the rule. Today only ``PlotModel``
handles ``ClearStageEvent`` and ``StageSeriesEvent``, but staging is a concept
like focus, which the user directs and other models may come to depend on, so
it is broadcast.

Event catalog
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 22 30 24

   * - Event
     - Published by
     - Subscribed by
     - Purpose
   * - ``FocusEvent``
     - ``LoadRawScanPresenter``
     - ``TaviProjectModel``
     - Tree selection, as bare uuids
   * - ``ClearFocusEvent``
     - ``TaviProjectModel``
     - ``PlotModel``, ``PlotterPresenter``, ``FittingPresenter``,
       ``DataFilePresenter``
     - Drop the old selection, stage included
   * - ``FocusRawScanEvent``
     - ``TaviProjectModel``
     - ``PlotModel``, ``PlotterPresenter``
     - Add raw scans; ``PlotModel`` builds preview plots
   * - ``FocusPlotEvent``
     - ``TaviProjectModel``, ``PlotModel``
     - ``PlotModel``, ``PlotterPresenter``
     - Add plots to the focus and the canvas
   * - ``FocusFitEvent``
     - ``TaviProjectModel``
     - ``PlotModel``, ``PlotterPresenter``, ``FittingPresenter``
     - Add fits; their series are focused too
   * - ``SyncPlotEvent``
     - ``PlotModel``
     - ``PlotterPresenter``
     - Focused plots after an edit or removal
   * - ``ClearStageEvent``
     - ``PlotterPresenter``
     - ``PlotModel``
     - Empty the stage before restaging
   * - ``StageSeriesEvent``
     - ``PlotterPresenter``
     - ``PlotModel``
     - Add series to the stage
   * - ``SyncStageEvent``
     - ``PlotModel``
     - ``PlotterPresenter``, ``FittingPresenter``, ``DataFilePresenter``
     - Staged series and their scans; the first leads
   * - ``RecomputeFitEvent``
     - ``TaviProjectModel``
     - ``FitModel``
     - Recompute fits for display
   * - ``SyncFitEvent``
     - ``FitModel``
     - ``PlotterPresenter``, ``FittingPresenter``, ``FitWindowPresenter``
     - Computed fit curves and results
   * - ``SyncFitSpecEvent``
     - ``FitModel``
     - ``FittingPresenter``
     - One saved member's spec, without refitting
   * - ``SyncFitHistoryEvent``
     - ``TaviProjectModel``
     - ``FittingPresenter``
     - Whether a member can be undone or redone
   * - ``RestoreFitMemberEvent``
     - ``TaviProjectModel``
     - ``FitWindowPresenter``
     - A member was undone or redone
   * - ``SyncPeakParamsEvent``, ``SyncBackgroundParamsEvent``
     - ``FitModel``
     - ``FittingPresenter``
     - Suggested starting parameters
   * - ``SetFitComponentsVisibleEvent``
     - ``FittingPresenter``
     - ``PlotterPresenter``
     - Show or hide fit components
   * - ``SavePlotEvent``, ``SaveFitEvent``
     - ``PlotModel``, ``FitModel``
     - ``TaviProjectModel`` (and ``FitWindowPresenter`` for fits)
     - Write into ``TaviData``
   * - ``AddRawScanEvent``, ``AddPlotEvent``, ``AddFitEvent``
     - ``TaviProjectModel``
     - ``LoadRawScanPresenter``
     - Item entered the project
   * - ``RemoveRawScanEvent``, ``RemovePlotEvent``, ``RemoveFitEvent``
     - ``TaviProjectModel``
     - ``LoadRawScanPresenter``, and ``PlotModel``, ``FittingPresenter``,
       ``FitWindowPresenter`` where they hold that kind of item
     - Item left the project
   * - ``SyncRecentProjectsEvent``
     - ``TaviProjectModel``
     - ``FileMenuPresenter``
     - The recent-projects list
   * - ``ReportErrorEvent``
     - Models, ``WorkerPool``
     - ``RecoveryService``
     - Surface an error
   * - ``StartApplicationEvent``
     - ``MainPresenter``
     - ``TaviProjectModel``
     - Views are built; push initial state

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
