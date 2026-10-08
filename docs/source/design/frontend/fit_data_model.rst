Fit Data Model and Sequential Fitting
=====================================

Overview
--------

A fit, like a :doc:`Plot <plot_data_model>`, is a *collection*. It covers one or
more scans, and fitting one scan is just the case with one member. A
**sequential fit** is not a separate type. It is a fit whose members were fit in
order, each one starting from the previous member's fitted values. Once the fit
has run, nothing records that it was sequential: every member is recomputed from
its own saved spec, the same as any other fit.

This page covers the data model, how the *stage* decides between a sequential
fit and refitting one member, and the event flows involved. The stage is the
part of the focus that edits and fits apply to: every focused series while the
plotter's "Apply All" is checked, only the series picked in the "Current Plot"
dropdown while it isn't. ``PlotModel`` tracks it and syncs it out as
``SyncStageEvent``; the first staged series *leads*, and is the one the fitting
panel shows. The `Event map`_ shows every fit event and what it triggers.

Use cases covered
-----------------

* **Run a sequential fit.** Select several scans, set up the fitting panel for
  the first one (or leave values blank to have them guessed), leave "Apply All"
  checked so every focused series is staged, and press Perform Fit. Every staged
  series is fit in stage order, which is tree selection order. Each fit starts
  from the previous fit's result. Every member gets its own window, and the whole fit is
  saved in ``TaviData.fits`` as one ``FitEntry``.
* **Change one member's parameters.** Uncheck "Apply All", then pick a scan in
  the "Current Plot" dropdown, which stages it alone. The fitting panel reloads
  exactly the spec and result saved for that member. Edit it and press Perform Fit: only that member
  is refit, and the rest of the fit is left alone. TAVI does no validity check on
  any fit; judging whether a fit worked is up to the user.
* **Undo or redo a refit.** While dialing in one member, press Undo Fit to go
  back to that member's state (spec and result) before its last save, and Redo
  Fit to go forward again. The buttons act on the one staged series, and are
  enabled only while exactly one series is staged ("Apply All" unchecked, or
  only one series focused). History is kept per member for the session, not
  saved in the project.
* **Start over on a new selection.** Selecting anything in the project tree
  resets the panel's background and peak fields to their defaults, and the
  range follows the newly leading series. Undo/Redo disable until there is
  something to step to. Picking a raw scan again and fitting it starts a new
  fit with its own history. To carry on with an earlier fit and its history,
  select that fit in the tree.
* **Recover a saved fit.** Select the fit in the project tree. Every member is
  recomputed from its stored spec (not re-chained) and drawn on the main plotter
  only, over the data each member was fit against. No windows open, and open
  ones are not refreshed. Uncheck "Apply All" and use the dropdown to inspect
  each member.

Data model
----------

.. mermaid::

    classDiagram
        class FitSpec {
            +str range_min
            +str range_max
            +str background
            +ParamField background_constant
            +ParamField background_slope
            +list~PeakField~ peaks
        }
        class FitMember {
            +PlotSeries series
            +Optional~FitResultSummary~ result
            +source_scan_uuid()
        }
        class FitEntry {
            +UUID uuid
            +str name
            +list~FitMember~ members
            +member_for(source_scan_uuid)
            +with_members(members)
            +without_scans(scan_uuids)
        }
        FitSpec <|-- FitMember
        FitEntry "1" o-- "1..*" FitMember

* ``FitSpec`` is the fitting panel's raw field state with no series attached,
  so one spec can be applied to a whole batch.
* ``FitMember`` binds a spec to the one series it was fit against. It also keeps
  the fitted **parameters** (``result``), because the use case asks for them to
  be saved in the project. It never keeps the evaluated curve or the data: those
  are recomputed whenever they are shown.
* A member's identity inside its fit is its ``source_scan_uuid``, so it needs no
  uuid of its own. This is the same key the plotter already uses for curves.
* ``FitEntry.without_scans`` mirrors ``Plot.without_scans``: removing a scan
  drops only that member, and the fit is deleted only once no member is left.

Key decisions
-------------

**One type, not a SequentialFit beside FitEntry.** A separate type would have
needed its own events, tree node, purge rule and recompute path. All of those
would be copies of the single-fit ones, and after the run the two kinds of fit
behave the same. Seeding is a property of *running* the fit
(``FitRequest.seed_from_previous``), so it lives on the request, not on the
saved object.

**Focus, stage, calculate, sync.** Fit events follow the same pattern as the
rest of the widgets. A *focus* event says what the user is looking at
(``FocusFitEvent``), the *stage* says what of that the next fit applies to
(``SyncStageEvent``), a model *calculates* (``FitModel``), and a *sync* event
brings the UI in line with the result: ``SyncFitEvent`` for computed curves,
``SyncFitSpecEvent`` for a saved member's spec. Sync events are display only.
Writing to the project is a separate *save* event, ``SaveFitEvent``, mirroring
``SavePlotEvent``.

**Focusing a fit focuses its data.** ``FocusFitEvent`` only adds to the focus,
like every ``Focus*`` event, since ``ClearFocusEvent`` has already cleared the
old selection. ``PlotModel`` answers it by publishing a ``FocusPlotEvent`` for
each member's series that isn't already focused, so a fit's data is drawn
under its curve whether the fit was selected alone or alongside its scan, and
the plotter stages those series like any others. ``PlotterPresenter`` and
``FittingPresenter`` only mark the fit pending, so its recomputed curve and
result are shown when they arrive.

**The fitting panel works on the stage.** ``FittingPresenter`` keeps the staged
series from ``SyncStageEvent`` and fits exactly those, in stage order: several
staged series run a sequential fit, one is fit alone. It no longer tracks the
focus or the "Apply All" checkbox itself; staging already encodes both.

**Save and sync are separate events.** Perform Fit publishes ``SaveFitEvent``
and then ``SyncFitEvent``. Recomputing a saved fit for display publishes only
``SyncFitEvent``, because nothing about the fit changed. One combined event
would need a flag saying which of the two it was, and every subscriber that
cares (the project model, the fit windows) would have to branch on it. With two
events, each subscriber registers for the one it needs.

**Save and sync name a fit and carry some of its members.** Refitting one member
publishes only that member. ``TaviProjectModel`` merges each saved member into
the stored entry by scan (``FitEntry.with_members``), so the members an event
leaves out stay unchanged rather than being removed.

**The model resolves the data.** ``FitRequest`` carries series pointers, not
x/y/err arrays. ``FitModel`` resolves each one against its raw scan handle, as it
already did for recompute. A presenter has no way to hold data for every scan in
a batch, and shouldn't.

**Sync the spec instead of recomputing when the lead changes.** When a
``SyncStageEvent`` brings a new series to the lead and the panel knows a fit
covering it, ``FittingPresenter`` calls
``FitModel.sync_fit_spec(fit_uuid, scan_uuid)`` through the model proxy.
``FitModel`` reads the saved member from its live handle on ``TaviData.fits``
(the same way ``PlotModel`` holds the plots handle) and publishes
``SyncFitSpecEvent``, so nothing is refit. This is a direct model call, not a
request event, because only ``FitModel`` answers it. A field edit resyncs the
stage with the same lead, which changes nothing here. Inside a focus chain the
scan-to-fit map has just been cleared, so this never fires there and adds
nothing to the broker's depth budget (see :doc:`../../guides/event_broker`).

**Undo/redo is per member, and a step is one save.** Every change to a fit
arrives as a ``SaveFitEvent``, so ``TaviProjectModel`` records the member each
merge displaces in a ``FitHistory`` keyed by ``(fit uuid, source scan uuid)``.
Field edits that were never fit are not steps: the panel's line edits have
their own undo for typing, and the use case is going back to a *fit*. A step
from a sequential run is recorded per member too, so a bad run can be undone
one scan at a time, but only with that one series staged - undo is for dialing
one member in. Other choices:

* History is linear. Saving after an undo clears that member's redo stack.
* Each member keeps at most ``FIT_HISTORY_LIMIT`` (10) states, and the oldest
  is dropped first.
* A failed fit (``result=None``) is not recorded as a step, since there is no
  fit to go back to.
* History lives beside ``TaviData`` rather than in it, so dialing a fit in
  doesn't grow the saved project. It is forgotten when its fit or scan is
  removed.
* Undo and redo are direct calls, ``TaviProjectInterface.undo_fit_member`` and
  ``redo_fit_member``, made through the project proxy rather than published as
  events: only ``TaviProjectModel`` holds the history, so there is nothing to
  broadcast until it has changed something. What changed is then broadcast as
  ``SyncFitHistoryEvent``, ``RestoreFitMemberEvent`` and ``RecomputeFitEvent``.
* A restored member is recomputed from its spec through ``RecomputeFitEvent``,
  never replayed from a cached curve.
* ``FittingPresenter`` keeps only the ``(can undo, can redo)`` flags from
  ``SyncFitHistoryEvent`` for each member, never the states themselves.

**Each selection starts over.** ``TaviProjectModel`` publishes
``ClearFocusEvent`` first thing on every ``FocusEvent``, before the chain
that event sets off. A new selection is two operations, a clear and a focus,
so the clear is its own event and says nothing about any widget; each
subscriber decides what letting go of the old selection means for it.
``FittingPresenter`` resets the panel's fields and forgets
which fit covers each scan, so the panel never shows one scan's parameters
against another scan's data. Because the scan-to-fit map is cleared, Perform
Fit on a raw scan picked again mints a new fit with its own history. A fit
picked from the tree is mapped back by the ``SyncFitEvent`` that follows, and
its history (still held by ``TaviProjectModel``) comes back with it.
The stage is part of the focus, so it is cleared with it and the panel's
staged series are forgotten too. ``FocusPlotEvent`` can't be the cue, because
a mixed selection publishes several of them, and a ``FittingPresenter`` handler
on ``FocusEvent`` itself would run after the chain or before it depending on
registration order. ``SyncPlotEvent``, which ``PlotModel`` publishes for every
field edit or removal, is a refresh of the same focus and never resets the
panel.

**Windows only for a freshly saved fit.** ``FitWindowPresenter`` only touches
windows for a fit it has just seen a ``SaveFitEvent`` for. On the save, it notes
which windows the following ``SyncFitEvent`` should fill: every member of a
multi-member fit (opening any that are missing), or only the windows already
open for a one-member refit. The sync then consumes that note. A recompute is
synced but never saved, so browsing a fit from the project tree shows it on the
main plotter only and leaves open windows as they were. Without that rule, every
reselection would spawn a window per scan, and refitting one member would reopen
a window the user closed. Each window has a "Close All" button that closes every
fit window. Open windows are tracked by ``(fit uuid, source scan uuid)`` and
closed when their fit or scan leaves the project.

**Failures stay in the fit.** A member whose fit raised (bad range, unsupported
shape) is kept with ``result=None``, so the user can fix just that scan, and
seeding continues from the last member that succeeded. When several series are
fit together, their validation errors are reported as one message, not one
dialog per scan. If no member succeeds, nothing is saved.

Event map
---------

Every event that takes part in fitting, who publishes it, and what each
subscriber does with it. Rounded boxes are user actions, rectangles are events,
and each arrow out of an event is one subscriber. An event a handler publishes
in turn is drawn as a further arrow, so following any path from a user action
gives the full chain it sets off. Dotted arrows are follow-ups that only fire
under the stated condition, or that end in UI state rather than another event.

.. mermaid::

    flowchart LR
        classDef user fill:#fff3cd,stroke:#b8860b
        classDef event fill:#dbe9ff,stroke:#3a6ea5

        tree([Select items in project tree]):::user
        combo([Pick series in Current Plot]):::user
        applyall([Toggle Apply All]):::user
        perform([Perform Fit]):::user
        undoredo([Undo / Redo Fit]):::user
        separately([Toggle Plot Separately]):::user
        remove([Remove items from tree]):::user

        FocusEvent:::event
        ClearFocusEvent:::event
        FocusRawScanEvent:::event
        FocusPlotEvent:::event
        FocusFitEvent:::event
        RecomputeFitEvent:::event
        ClearStageEvent:::event
        StageSeriesEvent:::event
        SyncStageEvent:::event
        SaveFitEvent:::event
        SyncFitEvent:::event
        SyncFitSpecEvent:::event
        AddFitEvent:::event
        SetFitComponentsVisibleEvent:::event
        RemoveRawScanEvent:::event
        RemoveFitEvent:::event
        SyncFitHistoryEvent:::event
        RestoreFitMemberEvent:::event

        tree -- LoadRawScanPresenter --> FocusEvent
        FocusEvent -- "TaviProjectModel: always, before routing" --> ClearFocusEvent
        ClearFocusEvent -. "FittingPresenter: reset panel, forget scan-to-fit map and stage" .-> panel
        ClearFocusEvent -. "PlotModel: forget focus and stage" .-> lastplots[(focus and stage)]
        FocusEvent -- "TaviProjectModel: scans selected" --> FocusRawScanEvent
        FocusEvent -- "TaviProjectModel: plots selected" --> FocusPlotEvent
        FocusEvent -- "TaviProjectModel: fits selected, or attached to a selected plot (1st)" --> FocusFitEvent
        FocusEvent -- "TaviProjectModel: same fits (2nd)" --> RecomputeFitEvent
        FocusRawScanEvent -- "PlotModel: build previews" --> FocusPlotEvent
        FocusFitEvent -- "PlotModel: member series not yet focused" --> FocusPlotEvent
        FocusFitEvent -. "PlotterPresenter + FittingPresenter: mark fits pending" .-> panel
        FocusPlotEvent -. "PlotModel: add to focus" .-> lastplots
        FocusPlotEvent -- "PlotterPresenter: draw, then stage new series" --> StageSeriesEvent
        RecomputeFitEvent -- "FitModel: recompute each member" --> SyncFitEvent

        combo -- "PlotterPresenter: Apply All off (1st)" --> ClearStageEvent
        combo -- "PlotterPresenter: picked series (2nd)" --> StageSeriesEvent
        applyall -- "PlotterPresenter (1st)" --> ClearStageEvent
        applyall -- "PlotterPresenter: all focused, or the lead (2nd)" --> StageSeriesEvent
        ClearStageEvent -. "PlotModel: empty the stage" .-> lastplots
        StageSeriesEvent -- "PlotModel: add to stage" --> SyncStageEvent
        SyncStageEvent -. "FittingPresenter: Perform Fit covers these; new lead with a fit calls FitModel.sync_fit_spec" .-> syncspec[FitModel.sync_fit_spec]
        syncspec -- "member found" --> SyncFitSpecEvent
        SyncFitSpecEvent -. "FittingPresenter: still the lead series" .-> panel

        perform -- "FittingPresenter calls FitModel.perform_fit" --> performfit[FitModel.perform_fit]
        performfit -- "at least one member succeeded (1st)" --> SaveFitEvent
        performfit -- "same members (2nd)" --> SyncFitEvent

        SaveFitEvent -. "TaviProjectModel: merge members; new uuid only" .-> AddFitEvent
        SaveFitEvent -. "FitWindowPresenter: note windows to fill" .-> windows[(fit windows)]
        SaveFitEvent -- "TaviProjectModel: record displaced members" --> SyncFitHistoryEvent
        SyncFitHistoryEvent -. "FittingPresenter: enable Undo/Redo for the one staged member" .-> panel

        undoredo -- "FittingPresenter calls TaviProjectModel.undo_fit_member / redo_fit_member" --> undocall[TaviProjectModel.undo_fit_member]
        undocall -- "swap member (1st)" --> SyncFitHistoryEvent
        undocall -- "(2nd)" --> RestoreFitMemberEvent
        undocall -- "restored member only (3rd)" --> RecomputeFitEvent
        RestoreFitMemberEvent -. "FitWindowPresenter: note window, only if open" .-> windows
        AddFitEvent -- LoadRawScanPresenter --> treenode[(tree /Fits node)]
        SyncFitEvent -. "PlotterPresenter: append curves" .-> canvas[(main plotter)]
        SyncFitEvent -. "FittingPresenter: map staged scans to fit, load lead's member" .-> panel[(fitting panel)]
        SyncFitEvent -. "FitWindowPresenter: noted windows only" .-> windows

        separately -- FittingPresenter --> SetFitComponentsVisibleEvent
        SetFitComponentsVisibleEvent -. "PlotterPresenter: show/hide components" .-> canvas

        remove -- "TaviProjectModel.remove_items, TaviData.purge: scan removed" --> RemoveRawScanEvent
        remove -- "fit removed, or its last scan" --> RemoveFitEvent
        RemoveRawScanEvent -. "FitWindowPresenter: close that scan's windows" .-> windows
        RemoveFitEvent -. "FitWindowPresenter: close fit's windows; FittingPresenter: forget uuid; tree: drop node" .-> windows

Rectangles named ``FitModel.*`` or ``TaviProjectModel.*`` are direct model
calls made through the model proxy, not events. They run on the proxy's worker
thread, and whatever they publish is dispatched from there.

Nesting depth matters (the broker allows 5, see
:doc:`../../guides/event_broker`). The deepest fit chain is browsing a fit,
which reaches the limit: ``FocusEvent`` (1), ``FocusFitEvent`` (2),
``FocusPlotEvent`` from ``PlotModel`` (3), ``StageSeriesEvent`` from
``PlotterPresenter`` (4), and ``SyncStageEvent`` (5). Selecting raw scans runs
the same chain through ``FocusRawScanEvent``. No ``SyncStageEvent`` handler may
publish an event of its own: ``FittingPresenter``'s call to
``FitModel.sync_fit_spec`` is a proxy call, and inside a focus chain it never
happens, because the scan-to-fit map has just been cleared. The recompute runs
``RecomputeFitEvent`` (2) and ``SyncFitEvent`` (3), and never saves, so no
``AddFitEvent`` follows. Perform Fit, ``sync_fit_spec`` and undo/redo start from
a model call rather than an event, so ``SaveFitEvent``, ``SyncFitEvent``,
``SyncFitSpecEvent``, ``SyncFitHistoryEvent``, ``RestoreFitMemberEvent`` and
the undo's ``RecomputeFitEvent`` are at depth 1, ``AddFitEvent`` and the undo's
``SyncFitEvent`` at 2. Restaging from the dropdown or "Apply All" runs
``ClearStageEvent`` / ``StageSeriesEvent`` (1) and ``SyncStageEvent`` (2).
Adding a new link to any of these chains should be checked against that
budget.

Orderings relied on rather than enforced by the broker:

* ``TaviProjectModel`` publishes ``ClearFocusEvent`` *before* routing a
  ``FocusEvent``, so every subscriber has let go of the old selection (and its
  stage) before the new selection's ``Focus*``, ``SyncStageEvent`` and
  ``SyncFitEvent`` reach it.
* ``TaviProjectModel`` publishes ``FocusFitEvent`` *before*
  ``RecomputeFitEvent``, so the presenters have marked the fit pending by the
  time its ``SyncFitEvent`` arrives, and the fit's series are already focused
  and staged.
* ``FitModel.perform_fit`` publishes ``SaveFitEvent`` *before*
  ``SyncFitEvent``. ``FitWindowPresenter`` only fills windows it noted on the
  save, so this order is what tells it the sync is a fresh fit.
* ``TaviProjectModel`` publishes ``RestoreFitMemberEvent`` *before* the
  ``RecomputeFitEvent`` for an undo/redo, for the same reason: it is how
  ``FitWindowPresenter`` knows that one recompute should refresh an open window.

One ordering is deliberately *not* relied on: ``PlotterPresenter`` stages from
its own ``FocusPlotEvent`` handler, which may run before or after
``PlotModel``'s. ``PlotModel`` keeps a stage request for a series it hasn't
focused yet and syncs the stage once that series' ``FocusPlotEvent`` reaches it,
so either order converges on the same ``SyncStageEvent``.

Run a sequential fit
--------------------

.. mermaid::

    sequenceDiagram
        participant User
        participant FittingView
        participant FittingPresenter
        participant FitModel
        participant EventBroker
        participant TaviProjectModel
        participant PlotterPresenter
        participant FitWindowPresenter
        participant LoadRawScanPresenter

        Note over FittingPresenter: SyncStageEvent staged s1..sn ("Apply All" checked)
        User ->> FittingView: Perform Fit
        FittingView ->> FittingPresenter: perform_fit_clicked
        FittingPresenter ->> FitModel: perform_fit(FitRequest(spec, series=[s1..sn], seed_from_previous=True))
        loop each series, in stage order
            FitModel ->> FitModel: resolve x/y/err, fit, seed next spec from this result
        end
        FitModel ->> EventBroker: SaveFitEvent(fit_uuid, members=[m1..mn])
        EventBroker ->> TaviProjectModel: store FitEntry
        TaviProjectModel ->> EventBroker: AddFitEvent (new uuid only)
        EventBroker ->> LoadRawScanPresenter: add node under /Fits
        EventBroker ->> FitWindowPresenter: note a window per member
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes=[o1..on])
        EventBroker ->> PlotterPresenter: draw each curve on the main canvas
        EventBroker ->> FittingPresenter: remember fit_uuid per staged series, show the lead member
        EventBroker ->> FitWindowPresenter: open the noted windows

Browse a saved fit
------------------

Selecting a fit in the project tree shows it, over the data each member was fit
against, on the main plotter only. The recompute publishes ``SyncFitEvent`` with
no ``SaveFitEvent`` before it, so ``FitWindowPresenter`` has nothing noted: no
window opens and any open window is left as it was.

.. mermaid::

    sequenceDiagram
        participant User
        participant LoadRawScanPresenter
        participant EventBroker
        participant TaviProjectModel
        participant PlotterPresenter
        participant FittingPresenter
        participant PlotModel
        participant FitModel
        participant FitWindowPresenter

        User ->> LoadRawScanPresenter: select fit in tree
        LoadRawScanPresenter ->> EventBroker: FocusEvent(ids)
        EventBroker ->> TaviProjectModel: resolve ids to FitEntry
        TaviProjectModel ->> EventBroker: ClearFocusEvent
        EventBroker ->> PlotModel: forget focus and stage
        EventBroker ->> PlotterPresenter: empty canvas and dropdown
        EventBroker ->> FittingPresenter: reset panel, forget scan-to-fit map
        TaviProjectModel ->> EventBroker: FocusFitEvent(fits, scans)
        EventBroker ->> PlotterPresenter: mark fit pending
        EventBroker ->> FittingPresenter: mark fit pending
        EventBroker ->> PlotModel: member series not yet focused
        PlotModel ->> EventBroker: FocusPlotEvent(one plot per member series)
        EventBroker ->> PlotModel: add to focus
        EventBroker ->> PlotterPresenter: draw the data, list each series
        PlotterPresenter ->> EventBroker: StageSeriesEvent(all, or the first)
        EventBroker ->> PlotModel: add to stage
        PlotModel ->> EventBroker: SyncStageEvent(staged series, scans)
        EventBroker ->> FittingPresenter: Perform Fit now covers these
        TaviProjectModel ->> EventBroker: RecomputeFitEvent(fits)
        EventBroker ->> FitModel: recompute every member from its own spec
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes)
        EventBroker ->> PlotterPresenter: append curves
        EventBroker ->> FittingPresenter: map members to fit, show the lead member
        EventBroker -->> FitWindowPresenter: ignored (nothing noted, no save)

Selecting a saved plot follows the same path for every fit stamped on it
(``Plot.fits``). The plot's own ``FocusPlotEvent`` has already focused its
series, so ``PlotModel`` finds them focused and adds nothing for the fit, and
the curves overlay the plot.

Refit one member
----------------

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant FittingPresenter
        participant EventBroker
        participant PlotModel
        participant TaviProjectModel
        participant FitModel
        participant FitWindowPresenter

        User ->> Plot1DView: unchecks "Apply All"
        PlotterPresenter ->> EventBroker: ClearStageEvent, StageSeriesEvent(lead)
        EventBroker ->> PlotModel: restage
        PlotModel ->> EventBroker: SyncStageEvent([lead])
        EventBroker ->> FittingPresenter: Perform Fit now covers the lead alone
        User ->> Plot1DView: picks scan k in Current Plot
        PlotterPresenter ->> EventBroker: ClearStageEvent, StageSeriesEvent([scan k])
        EventBroker ->> PlotModel: restage
        PlotModel ->> EventBroker: SyncStageEvent([series k])
        EventBroker ->> PlotterPresenter: sync axis fields, point dropdown at scan k
        EventBroker ->> FittingPresenter: scan k leads, and a fit is known for it
        FittingPresenter ->> FitModel: sync_fit_spec(fit_uuid, scan k)
        FitModel ->> EventBroker: SyncFitSpecEvent(fit_uuid, member k)
        EventBroker ->> FittingPresenter: load member k's spec + result into the panel
        User ->> FittingPresenter: edits, Perform Fit
        FittingPresenter ->> FitModel: perform_fit(FitRequest(series=[scan k], fit_uuid))
        FitModel ->> EventBroker: SaveFitEvent(fit_uuid, members=[mk])
        EventBroker ->> TaviProjectModel: merge member k, others untouched
        EventBroker ->> FitWindowPresenter: note scan k's window, only if open
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes=[ok])
        EventBroker ->> PlotterPresenter: replace scan k's curve
        EventBroker ->> FitWindowPresenter: refresh it if noted (never reopens)

Undo a refit
------------

.. mermaid::

    sequenceDiagram
        participant User
        participant FittingView
        participant FittingPresenter
        participant EventBroker
        participant TaviProjectModel
        participant FitWindowPresenter
        participant FitModel
        participant PlotterPresenter

        Note over TaviProjectModel: each SaveFitEvent recorded the member it displaced
        User ->> FittingView: Undo Fit (enabled: one series staged, and it has a step)
        FittingView ->> FittingPresenter: undo_fit_clicked
        FittingPresenter ->> TaviProjectModel: undo_fit_member(fit_uuid, scan k) via proxy
        TaviProjectModel ->> TaviProjectModel: pop member k's undo state, push current onto redo
        TaviProjectModel ->> TaviProjectModel: fits[fit_uuid].with_members([restored])
        TaviProjectModel ->> EventBroker: SyncFitHistoryEvent(can_undo, can_redo)
        EventBroker ->> FittingPresenter: update Undo/Redo buttons
        TaviProjectModel ->> EventBroker: RestoreFitMemberEvent(fit_uuid, scan k)
        EventBroker ->> FitWindowPresenter: note scan k's window, only if open
        TaviProjectModel ->> EventBroker: RecomputeFitEvent(fit with member k only)
        EventBroker ->> FitModel: recompute member k from its restored spec
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes=[ok])
        EventBroker ->> PlotterPresenter: replace scan k's curve
        EventBroker ->> FittingPresenter: load restored member into the panel
        EventBroker ->> FitWindowPresenter: refresh it if noted

Redo is the same with the stacks swapped.

Remove a fit or scan
--------------------

.. mermaid::

    sequenceDiagram
        participant User
        participant LoadRawScanPresenter
        participant TaviProjectModel
        participant EventBroker
        participant FittingPresenter
        participant FitWindowPresenter

        User ->> LoadRawScanPresenter: remove items
        LoadRawScanPresenter ->> TaviProjectModel: remove_items(uuids)
        TaviProjectModel ->> TaviProjectModel: TaviData.purge - drop members of removed scans, fits left empty
        TaviProjectModel ->> TaviProjectModel: forget undo history of removed fits and scans
        TaviProjectModel ->> EventBroker: RemoveRawScanEvent per scan
        EventBroker ->> FitWindowPresenter: close every window on that scan
        EventBroker ->> LoadRawScanPresenter: drop tree node
        TaviProjectModel ->> EventBroker: RemoveFitEvent per fit removed or emptied
        EventBroker ->> FittingPresenter: forget fit uuid (next Perform Fit mints a new one)
        EventBroker ->> FitWindowPresenter: close the fit's windows
        EventBroker ->> LoadRawScanPresenter: drop tree node
