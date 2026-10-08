Selection → Visualization Flow
==============================

Overview
--------

Selecting items in the project tree triggers a chain of typed events that
ends with the selection drawn in the ``Plot1DView`` widget. The chain is fully
decoupled: no layer holds a direct reference to any other layer. Each
component publishes an event and relies on the ``EventBroker`` singleton to
deliver it to the next stage.

Two concepts run through every flow on this page:

- **Focus** is what is being looked at: every series drawn on the canvas and
  listed in the "Current Plot" dropdown.
- **Stage** is the part of the focus that edits and fits apply to. With
  "Apply All" checked every focused series is staged; without it, only the
  series picked in the "Current Plot" dropdown is. The first staged series
  *leads*: it is the one the plotter's axis/preset fields, the Data File tab
  and the fitting panel show.

Both are changed by small, single-purpose events. Each one does exactly one
thing, and a compound action is built from them:

- A new tree selection is ``ClearFocusEvent`` followed by one or more
  ``Focus*`` events. ``Focus*`` events only *add* to what is focused, so a
  mixed scan + plot + fit selection simply publishes several of them side by
  side.
- A restage is ``ClearStageEvent`` followed by ``StageSeriesEvent``, which
  likewise only adds.

The tree supports multi-select: every selected item is focused, in the order
the user selected them (see `Selection order is preserved across
multi-select`_ below).


Participants
------------

- **Project tree** — translates user interaction into a generic selection
  event, preserving the order items were selected in
- **Load presenter** — collects the selected identifiers and publishes
  ``FocusEvent``
- **Project model** (``TaviProjectModel``) — the only component that knows
  what type each tree uuid is. It clears the old focus, resolves each
  identifier to a domain object, and publishes one additive ``Focus*`` event
  per item type
- **Plot model** (``PlotModel``) — owns the focus (``_last_plots``) and the
  stage (``_staged_uuids``). It turns raw scans and fits into single-series
  preview ``Plot``\ s, applies field edits to the staged series, and syncs both
  focus and stage back out (``SyncPlotEvent``, ``SyncStageEvent``)
- **Plotter presenter + view** — draws focused series, lists them in the
  "Current Plot" dropdown, and publishes the user's staging choices
  (``ClearStageEvent``/``StageSeriesEvent``)
- **Data file presenter + view** — shows the lead staged series' scan in the
  Data File tab
- **Fitting presenter** — fits the staged series and shows the lead one in
  the fitting panel (see :doc:`fit_data_model`)

Neither model ever hands a live ``raw_scans``/``plots`` handle to a
presenter. Whatever a presenter or view needs to render travels *inside* the
event that triggers it — see :doc:`plot_data_model` for why.


Event summary
-------------

=========================  =========================  ==========================================================
Event                      Published by               Meaning
=========================  =========================  ==========================================================
``FocusEvent``             ``LoadRawScanPresenter``   The tree selection changed; ids of every selected item.
``ClearFocusEvent``        ``TaviProjectModel``       Nothing is focused any more (stage included).
``FocusRawScanEvent``      ``TaviProjectModel``       Add these raw scans to the focus.
``FocusPlotEvent``         ``TaviProjectModel``,      Add these plots to the focus.
                           ``PlotModel``
``FocusFitEvent``          ``TaviProjectModel``       Add these fits (and the series they were fit against).
``SyncPlotEvent``          ``PlotModel``              The focused plots as they now stand, after an edit or
                                                      removal. Same focus, new content.
``ClearStageEvent``        ``PlotterPresenter``       Nothing is staged any more.
``StageSeriesEvent``       ``PlotterPresenter``       Add these focused series to the stage.
``SyncStageEvent``         ``PlotModel``              The staged series and their scans; the first one leads.
=========================  =========================  ==========================================================


Selection — Full Chain
----------------------

The diagram shows a mixed selection of raw scans and saved plots with "Apply
All" checked. A selection of only one type simply skips the other branch.

.. mermaid::

    sequenceDiagram
        participant User
        participant LoadRawScanPresenter
        participant EventBroker
        participant TaviProjectModel
        participant PlotModel
        participant PlotterPresenter
        participant Plot1DView
        participant DataFilePresenter

        User ->> LoadRawScanPresenter: selects scan and plot nodes (in selection order)
        LoadRawScanPresenter ->> EventBroker: FocusEvent(ids)

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> EventBroker: ClearFocusEvent
        EventBroker ->> PlotModel: forget focused plots and stage
        EventBroker ->> PlotterPresenter: empty canvas and dropdown, reset controls
        EventBroker ->> DataFilePresenter: clear the data tab

        TaviProjectModel ->> TaviProjectModel: resolve ids → RawScan / Plot / FitEntry
        TaviProjectModel ->> EventBroker: FocusRawScanEvent(scans)
        EventBroker ->> PlotterPresenter: offer the first scan's columns as preset channels
        EventBroker ->> PlotModel: build one single-series preview Plot per scan
        PlotModel ->> EventBroker: FocusPlotEvent(previews, scans)
        EventBroker ->> PlotModel: append previews to _last_plots
        EventBroker ->> PlotterPresenter: resolve_series against event.scans, add to canvas and dropdown
        PlotterPresenter ->> Plot1DView: add_plots_signal(resolved) — draws without clearing
        PlotterPresenter ->> EventBroker: StageSeriesEvent(new series uuids)
        EventBroker ->> PlotModel: add to stage
        PlotModel ->> EventBroker: SyncStageEvent(staged series, scans)
        EventBroker ->> PlotterPresenter: sync axis/preset fields to the lead series
        EventBroker ->> DataFilePresenter: show the lead series' scan

        TaviProjectModel ->> EventBroker: FocusPlotEvent(saved plots, scans)
        Note over EventBroker,DataFilePresenter: same handling as the preview FocusPlotEvent above

With "Apply All" unchecked, ``PlotterPresenter`` stages only the first series
of the selection — so there is always something for edits and fits to act
on — and leaves later ``FocusPlotEvent``\ s unstaged.

Selecting fits adds one more branch. ``TaviProjectModel`` publishes
``FocusFitEvent``, and ``PlotModel`` answers with a ``FocusPlotEvent`` for each
fitted series that is not already focused, so a fit's data is drawn under its
curve. The fit curves themselves arrive later through ``RecomputeFitEvent`` and
``SyncFitEvent`` — see :doc:`fit_data_model`.


Restaging — the dropdown and "Apply All"
----------------------------------------

Changing what is staged never re-renders the canvas: nothing that is drawn
has changed, only which series edits apply to.

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant EventBroker
        participant PlotModel
        participant DataFilePresenter
        participant FittingPresenter

        alt picks a series in "Current Plot" ("Apply All" off)
            User ->> Plot1DView: picks entry k
            Plot1DView ->> PlotterPresenter: plot_combo_index_changed(k)
            PlotterPresenter ->> EventBroker: ClearStageEvent
            PlotterPresenter ->> EventBroker: StageSeriesEvent([series k])
        else toggles "Apply All"
            User ->> Plot1DView: toggles the checkbox
            Plot1DView ->> PlotterPresenter: apply_all_toggled(checked)
            PlotterPresenter ->> EventBroker: ClearStageEvent
            PlotterPresenter ->> EventBroker: StageSeriesEvent(lead first, then every other focused series — or the lead alone)
        end
        EventBroker ->> PlotModel: replace stage
        PlotModel ->> EventBroker: SyncStageEvent(staged series, scans)
        EventBroker ->> PlotterPresenter: sync fields, point the dropdown at the lead
        EventBroker ->> DataFilePresenter: show the lead series' scan
        EventBroker ->> FittingPresenter: track the stage, show the lead's saved fit if it has one

The dropdown is disabled while "Apply All" is checked, so a dropdown pick
always stages one series alone. Checking "Apply All" keeps the current lead
first, so the panels keep showing the series they showed before.


Editing the staged series
-------------------------

Editing an axis or preset field applies to the staged series only. The
presenter does not pick a target: ``PlotModel`` already tracks the stage.

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant PlotModel
        participant EventBroker

        User ->> Plot1DView: edits an axis/preset field
        Plot1DView ->> PlotterPresenter: fields_focus_changed
        PlotterPresenter ->> PlotModel: update_fields(fields) (proxy call)
        PlotModel ->> PlotModel: update every staged series; unstaged series carried through
        PlotModel ->> EventBroker: SyncPlotEvent(all focused plots, scans)
        EventBroker ->> PlotterPresenter: redraw the canvas, re-append fit curves
        PlotModel ->> EventBroker: SyncStageEvent(staged series, scans)
        EventBroker ->> PlotterPresenter: sync fields to the edited lead

If any staged series rejects the fields (an unknown column, a non-numeric
preset value), ``PlotModel`` reports the error and changes nothing.

Removing a focused scan from the project takes the same path:
``PlotModel`` drops the scan's series and publishes ``SyncPlotEvent``, and
``SyncStageEvent`` too if a staged series went with it. If every staged
series was removed, ``PlotterPresenter`` restages the first remaining one.

The plotter's "Show Title" checkbox is a direct proxy call,
``PlotModel.set_show_title(show_title)``, not an event: only ``PlotModel``
owns series labels, so there is no one else to tell. It applies to the next
scans focused.


Add Plot — Saving a New Plot Entry
----------------------------------

Clicking **Add Plot** captures every focused series as one new, independent
plot in the project.

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant PlotModel
        participant EventBroker
        participant TaviProjectModel
        participant LoadRawScanPresenter

        User ->> Plot1DView: clicks "Add Plot"
        Plot1DView ->> PlotterPresenter: plot_clicked
        PlotterPresenter ->> PlotModel: save_focused_plots(fit_uuids of drawn curves)
        PlotModel ->> PlotModel: copy every focused series into one new Plot (fresh uuid)
        PlotModel ->> EventBroker: SavePlotEvent(plot)
        EventBroker ->> TaviProjectModel: store in TaviData.plots
        TaviProjectModel ->> EventBroker: AddPlotEvent(uuid, friendly_name, friendly_path="")
        EventBroker ->> LoadRawScanPresenter: add the node under /Plots

The new ``Plot`` always gets a fresh uuid, so clicking **Add Plot** twice on
the same focus creates two independent entries rather than overwriting one.
The fits currently drawn on the canvas are stamped onto it (``Plot.fits``),
so focusing it later brings its fit curves back too.


Key Design Decisions
--------------------

Prime events: clear, then add
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A new selection is two operations, a clear and a focus, so it is published as
two kinds of event rather than one "change selection" event. Every ``Focus*``
event only adds, which is what lets a mixed selection publish
``FocusRawScanEvent``, ``FocusPlotEvent`` and ``FocusFitEvent`` side by side
without one wiping the canvas the other just drew. Earlier versions folded
the clear into the focus events instead — a replacing ``PlotFocusEvent``, an
``exclusive`` flag on the fit event, and an ``also_plots`` field so a scan +
plot selection could be merged into a single publish. Splitting the clear out
removed all three.

The stage follows the same pattern (``ClearStageEvent`` then
``StageSeriesEvent``). ``ClearFocusEvent`` implies the stage is cleared too,
since the stage is part of the focus.

Focus versus sync
~~~~~~~~~~~~~~~~~

``FocusPlotEvent`` and ``SyncPlotEvent`` carry the same payload but mean
different things, and their subscribers react differently. A focus is a new
selection: it may reset controls and pick a default stage. A sync is the same
focus with updated content: it redraws but must not reset anything the user
set up. ``PlotModel`` publishes a sync after every field edit, so a single
event for both would wipe the panels on every keystroke.

The user stages; the model tracks and syncs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Staging is a user choice, made in the plotter's UI, so ``PlotterPresenter``
publishes it — the same way the tree publishes ``FocusEvent``. ``PlotModel``
is the source of truth for what is actually staged: it resolves the staged
uuids against the focused plots and publishes ``SyncStageEvent`` with the
series and scans. Every consumer reads the stage from that one sync, so the
data tab, the plotter fields and the fitting panel can never disagree about
which series leads.

``StageSeriesEvent`` may arrive before the plots it names: ``PlotterPresenter``
stages from its own ``FocusPlotEvent`` handler, which can run before
``PlotModel``'s. ``PlotModel`` keeps the uuids anyway and syncs the stage
again once a ``FocusPlotEvent`` brings the series in, so the result does not
depend on subscriber registration order.

There is no separate "active series". With "Apply All" off the stage is
exactly the series picked in the dropdown; with it on the dropdown is
disabled and the lead staged series plays that role.

Events for broadcast, proxy calls for requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Focus, stage and sync are broadcast as events because many components react
to each. A request that only one component answers — a field edit
(``update_fields``), "Show Title" (``set_show_title``), Add Plot
(``save_focused_plots``) — is a direct call through that model's proxy
instead. The model then publishes a sync event if the data changed.

Typed event routing in the project model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only ``TaviData`` knows what type each tree uuid is, so the tree publishes
one generic ``FocusEvent`` and ``TaviProjectModel`` narrows it into typed
``Focus*`` events. The tree never needs to know whether it selected a scan, a
plot or a fit.

PlotModel as an adapter
~~~~~~~~~~~~~~~~~~~~~~~

``PlotModel`` converts raw scans (and fitted series) into ``Plot``
compositions (``Plot``/``PlotSeries`` — see :doc:`plot_data_model`). A
multi-scan focus produces one single-series preview ``Plot`` per scan, not one
multi-series ``Plot``, so each run stays independently stageable in the
"Current Plot" dropdown. Preview plots are never written into
``TaviData.plots``; only **Add Plot** persists them.

Series are identified by their source scan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dropdown, the stage and the drawn fit curves all key a series by its
``source_scan_uuid``, never by its containing ``Plot``'s uuid. That is what
lets one series be picked out of a fused, multi-series saved plot exactly as
it would be among several single-series previews, and what lets a scan that
is selected alongside its own fit be focused once rather than twice.

Events carry their own data; presenters/views hold none
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``FocusPlotEvent``, ``SyncPlotEvent`` and ``SyncStageEvent`` carry the
``Scan`` objects their series reference (``scans``), gathered by the model
that publishes them. ``PlotterPresenter`` resolves each series against that
snapshot and forwards resolved arrays to the view. Neither the presenter nor
the view ever holds a live handle to ``raw_scans``; between events the
presenter keeps only uuids and dropdown labels.

Controls are reset on clear, synced from the lead
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **ClearFocusEvent** — ``PlotterPresenter`` calls
  ``reset_controls_to_defaults()``: rebin back to *No Rebin* with
  ``0``/``2``/``0.02``, preset type back to ``NONE``, preset value back to
  ``"1"``. It does **not** derive normalization from
  ``scan.tavimeta.normalization``.
- **FocusRawScanEvent** — the preset channel dropdown is repopulated from
  the first focused scan's columns.
- **SyncStageEvent** — the axis and preset fields are synced to the lead
  staged series (``sync_axis_fields``/``sync_preset_fields``), so they always
  describe the series an edit would change.

Both sync methods and the reset block widget signals while writing, so
programmatically updating a control never re-emits ``fields_focus_changed``
and loops back into the model.

Selection order is preserved across multi-select
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``TreeViewWidget`` tracks a ``_selection_order`` list alongside Qt's own
selection model, appending a uuid when it's selected and removing it when
it's deselected. ``get_selected_items()`` returns items in that order rather
than Qt's (unordered) ``selectedIndexes()`` order, and drops any uuid no
longer selected (e.g. a removed tree entry) before returning. Deselecting an
item and reselecting it moves it to the end. Within each item type this
determines the order series appear in the "Current Plot" dropdown and, with
"Apply All" off, which one is staged by default (the first). Across types,
raw scans are focused before saved plots, and fits last.

Depth budget
~~~~~~~~~~~~

The broker allows a nesting depth of 5 (see :doc:`../../guides/event_broker`).
The deepest chain here is a raw scan selection: ``FocusEvent`` (1),
``FocusRawScanEvent`` (2), ``FocusPlotEvent`` (3), ``StageSeriesEvent`` (4),
``SyncStageEvent`` (5). Handlers of ``SyncStageEvent`` reached this way must
not publish anything further.
