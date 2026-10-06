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

This page covers the data model, how the "Apply All" checkbox decides between a
sequential fit and refitting one member, and the event flows involved. The
`Event map`_ shows every fit event and what it triggers.

Use cases covered
-----------------

* **Run a sequential fit.** Select several scans, set up the fitting panel for
  the first one (or leave values blank to have them guessed), leave "Apply All"
  checked, and press Perform Fit. Every focused series is fit in "Current Plot"
  dropdown order, which is tree selection order. Each fit starts from the
  previous fit's result. Every member gets its own window, and the whole fit is
  saved in ``TaviData.fits`` as one ``FitEntry``.
* **Change one member's parameters.** Uncheck "Apply All", then pick a scan in
  the "Current Plot" dropdown. The fitting panel reloads exactly the spec and
  result saved for that member. Edit it and press Perform Fit: only that member
  is refit, and the rest of the fit is left alone. TAVI does no validity check on
  any fit; judging whether a fit worked is up to the user.
* **Recover a saved fit.** Select the fit in the project tree. Every member is
  recomputed from its stored spec (not re-chained) and drawn on the main plotter
  only. No windows open, and open ones are not refreshed. Uncheck "Apply All"
  and use the dropdown to inspect each member.

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

**Focus, calculate, sync.** Fit events follow the same pattern as the rest of
the widgets. A *focus* event says what the user selected (``FitFocusEvent``), a
model *calculates* (``FitModel``), and a *sync* event brings the UI in line with
the result: ``SyncFitEvent`` for computed curves, ``SyncFitSpecEvent`` for a
saved member's spec. Sync events are display only. Writing to the project is a
separate *save* event, ``SaveFitEvent``, mirroring ``SavePlotEvent``.

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

**Sync the spec instead of recomputing on a dropdown switch.** Picking a series
in the "Current Plot" dropdown makes ``FittingPresenter`` call
``FitModel.sync_fit_spec(fit_uuid, scan_uuid)`` through the model proxy.
``FitModel`` reads the saved member from its live handle on ``TaviData.fits``
(the same way ``PlotModel`` holds the plots handle) and publishes
``SyncFitSpecEvent``, so nothing is refit. This is a direct model call, not a
request event, because only ``FitModel`` answers it. ``FittingPresenter``
listens to the dropdown's own ``FocusActivePlotEvent`` rather than
``ActivePlotChangedEvent``. The latter fires deep inside every focus chain,
and syncing from there would eat into the broker's depth budget (see
:doc:`../../guides/event_broker`).

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
        separately([Toggle Plot Separately]):::user
        remove([Remove items from tree]):::user

        FocusEvent:::event
        RawScanFocusEvent:::event
        PlotFocusEvent:::event
        FitFocusEvent:::event
        FitRecomputeEvent:::event
        SaveFitEvent:::event
        SyncFitEvent:::event
        SyncFitSpecEvent:::event
        FitAppendEvent:::event
        ActivePlotChangedEvent:::event
        FocusActivePlotEvent:::event
        ApplyAllChangedEvent:::event
        FitComponentsVisibilityChangedEvent:::event
        RawScanRemoveEvent:::event
        FitRemoveEvent:::event

        tree -- LoadRawScanPresenter --> FocusEvent
        FocusEvent -- "TaviProjectModel: scans selected" --> RawScanFocusEvent
        FocusEvent -- "TaviProjectModel: plots only" --> PlotFocusEvent
        FocusEvent -- "TaviProjectModel: fits selected, or attached to a selected plot (1st)" --> FitFocusEvent
        FocusEvent -- "TaviProjectModel: same fits (2nd)" --> FitRecomputeEvent
        RawScanFocusEvent -- "PlotModel: build previews" --> PlotFocusEvent
        PlotFocusEvent -- "PlotterPresenter: render, reset dropdown" --> ActivePlotChangedEvent
        PlotFocusEvent -. "FittingPresenter: cache focused series" .-> focused[(focused series)]
        FitFocusEvent -- "PlotterPresenter: exclusive only, render fit's scans" --> ActivePlotChangedEvent
        FitFocusEvent -. "FittingPresenter: mark fits selected" .-> focused
        FitFocusEvent -. "PlotModel: sync _last_plots, no publish" .-> lastplots[(_last_plots)]
        FitRecomputeEvent -- "FitModel: recompute each member" --> SyncFitEvent

        perform -- "FittingPresenter calls FitModel.perform_fit" --> performfit[FitModel.perform_fit]
        performfit -- "at least one member succeeded (1st)" --> SaveFitEvent
        performfit -- "same members (2nd)" --> SyncFitEvent

        SaveFitEvent -. "TaviProjectModel: merge members; new uuid only" .-> FitAppendEvent
        SaveFitEvent -. "FitWindowPresenter: note windows to fill" .-> windows[(fit windows)]
        FitAppendEvent -- LoadRawScanPresenter --> treenode[(tree /Fits node)]
        SyncFitEvent -. "PlotterPresenter: append curves" .-> canvas[(main plotter)]
        SyncFitEvent -. "FittingPresenter: map scan to fit, load panel" .-> panel[(fitting panel)]
        SyncFitEvent -. "FitWindowPresenter: noted windows only" .-> windows

        combo -- PlotterPresenter --> FocusActivePlotEvent
        FocusActivePlotEvent -- "TaviProjectModel or PlotModel: owner resolves series" --> ActivePlotChangedEvent
        FocusActivePlotEvent -. "FittingPresenter: series has a fit, calls FitModel.sync_fit_spec" .-> syncspec[FitModel.sync_fit_spec]
        syncspec -- "member found" --> SyncFitSpecEvent
        SyncFitSpecEvent -. "FittingPresenter: still the active series" .-> panel
        ActivePlotChangedEvent -. "FittingPresenter + PlotterPresenter: track active series, sync fields" .-> panel

        applyall -- PlotterPresenter --> ApplyAllChangedEvent
        ApplyAllChangedEvent -. "FittingPresenter: Perform Fit covers all or one" .-> panel
        separately -- FittingPresenter --> FitComponentsVisibilityChangedEvent
        FitComponentsVisibilityChangedEvent -. "PlotterPresenter: show/hide components" .-> canvas

        remove -- "TaviProjectModel.remove_items, TaviData.purge: scan removed" --> RawScanRemoveEvent
        remove -- "fit removed, or its last scan" --> FitRemoveEvent
        RawScanRemoveEvent -. "FitWindowPresenter: close that scan's windows" .-> windows
        FitRemoveEvent -. "FitWindowPresenter: close fit's windows; FittingPresenter: forget uuid; tree: drop node" .-> windows

Rectangles named ``FitModel.*`` are direct model calls made through the model
proxy, not events. They run on the proxy's worker thread, and whatever they
publish is dispatched from there.

Nesting depth matters (the broker allows 5, see
:doc:`../../guides/event_broker`). The deepest fit chains are at depth 3:
browsing a fit runs ``FocusEvent`` (1), ``FitFocusEvent`` (2),
``ActivePlotChangedEvent`` (3), then ``FitRecomputeEvent`` (2),
``SyncFitEvent`` (3). A recompute never saves, so no ``FitAppendEvent``
follows. Perform Fit and ``sync_fit_spec`` start from a model call rather than
an event, so ``SaveFitEvent``, ``SyncFitEvent`` and ``SyncFitSpecEvent`` are at
depth 1 and ``FitAppendEvent`` at 2. Adding a new link to any of these chains
should be checked against that budget.

Two orderings are relied on rather than enforced by the broker:

* ``TaviProjectModel`` publishes ``FitFocusEvent`` *before*
  ``FitRecomputeEvent``, so the presenters have marked the fit selected by the
  time its ``SyncFitEvent`` arrives.
* ``FitModel.perform_fit`` publishes ``SaveFitEvent`` *before*
  ``SyncFitEvent``. ``FitWindowPresenter`` only fills windows it noted on the
  save, so this order is what tells it the sync is a fresh fit.

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

        User ->> FittingView: Perform Fit ("Apply All" checked, several series focused)
        FittingView ->> FittingPresenter: perform_fit_clicked
        FittingPresenter ->> FitModel: perform_fit(FitRequest(spec, series=[s1..sn], seed_from_previous=True))
        loop each series, in dropdown order
            FitModel ->> FitModel: resolve x/y/err, fit, seed next spec from this result
        end
        FitModel ->> EventBroker: SaveFitEvent(fit_uuid, members=[m1..mn])
        EventBroker ->> TaviProjectModel: store FitEntry
        TaviProjectModel ->> EventBroker: FitAppendEvent (new uuid only)
        EventBroker ->> LoadRawScanPresenter: add node under /Fits
        EventBroker ->> FitWindowPresenter: note a window per member
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes=[o1..on])
        EventBroker ->> PlotterPresenter: draw each curve on the main canvas
        EventBroker ->> FittingPresenter: remember fit_uuid per series, show the active member
        EventBroker ->> FitWindowPresenter: open the noted windows

Browse a saved fit
------------------

Selecting a fit in the project tree shows it on the main plotter only. The
recompute publishes ``SyncFitEvent`` with no ``SaveFitEvent`` before it, so
``FitWindowPresenter`` has nothing noted: no window opens and any open window is
left as it was.

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
        TaviProjectModel ->> EventBroker: FitFocusEvent(fits, exclusive, scans)
        EventBroker ->> PlotterPresenter: render each member's scan, reset dropdown, mark fit pending
        PlotterPresenter ->> EventBroker: ActivePlotChangedEvent(first series)
        EventBroker ->> FittingPresenter: track active series
        EventBroker ->> FittingPresenter: mark fit selected, cache its series for Apply All
        EventBroker ->> PlotModel: sync _last_plots (no publish)
        TaviProjectModel ->> EventBroker: FitRecomputeEvent(fits)
        EventBroker ->> FitModel: recompute every member from its own spec
        FitModel ->> EventBroker: SyncFitEvent(fit_uuid, outcomes)
        EventBroker ->> PlotterPresenter: append curves
        EventBroker ->> FittingPresenter: map series to fit, show active member
        EventBroker -->> FitWindowPresenter: ignored (nothing noted, no save)

Selecting a saved plot follows the same path for every fit stamped on it
(``Plot.fits``), except ``FitFocusEvent`` is non-exclusive: the plot's own
``PlotFocusEvent`` has already rendered the canvas, and the curves overlay it.

Refit one member
----------------

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant FittingPresenter
        participant EventBroker
        participant TaviProjectModel
        participant FitModel
        participant FitWindowPresenter

        User ->> Plot1DView: unchecks "Apply All"
        PlotterPresenter ->> EventBroker: ApplyAllChangedEvent(False)
        EventBroker ->> FittingPresenter: Perform Fit now covers the active series only
        User ->> Plot1DView: picks scan k in Current Plot
        PlotterPresenter ->> EventBroker: FocusActivePlotEvent(scan k)
        EventBroker ->> TaviProjectModel: resolve series (or PlotModel, whichever owns it)
        TaviProjectModel ->> EventBroker: ActivePlotChangedEvent(scan k)
        EventBroker ->> FittingPresenter: scan k is active
        EventBroker ->> FittingPresenter: fit known for scan k
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
        TaviProjectModel ->> EventBroker: RawScanRemoveEvent per scan
        EventBroker ->> FitWindowPresenter: close every window on that scan
        EventBroker ->> LoadRawScanPresenter: drop tree node
        TaviProjectModel ->> EventBroker: FitRemoveEvent per fit removed or emptied
        EventBroker ->> FittingPresenter: forget fit uuid (next Perform Fit mints a new one)
        EventBroker ->> FitWindowPresenter: close the fit's windows
        EventBroker ->> LoadRawScanPresenter: drop tree node


Open questions
--------------

* **Undo/redo.** Not implemented. Because every edit to a fit arrives as a
  ``SaveFitEvent`` naming the fit and the members it replaces, a natural
  design is a command stack in ``TaviProjectModel`` that records the members a
  merge displaced. Undo would re-merge those members and recompute them through
  ``FitRecomputeEvent``.
