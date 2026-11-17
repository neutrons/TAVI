Group Fit Window
================

Overview
--------

A **group fit** fits a sweep of scans (typically in temperature or field) in one
go, with :ref:`GroupFit <group_fit>` doing the fitting. It is separate from the
sequential fit in :doc:`fit_data_model`: it has its own window, its own model
and presenter, and it is display only. Nothing is saved into ``TaviData`` yet.

The Group Fit window lists the group's scans beside a canvas and a copy of the
1D fitting panel. Each scan keeps its own draft ``FitSpec`` in the window, and
one y-vs-x column pair is shared by the whole group.

Use cases covered
-----------------

* **Build a group.** Right-click scans or folders in the project tree and pick
  "Add to Group Fit". The window opens and the scans are appended in tree order.
  A scan already in the group keeps its place and spec.
* **Set the group up.** Edit the first scan in the panel. It is the group model,
  so each edit also reaches every other scan's field that still holds the old
  value. A field edited on another scan alone stays that scan's own change.
* **Guess starting values.** A peak's or the background's Suggest Params. guesses
  from the data of the scan on screen, over its own range, against the group's
  y-vs-x columns. A guess on the first scan reaches the others like any other
  edit; a guess on another scan is that scan's own change.
* **Fit and compare.** Press Group Fit, then select one scan or several
  (Ctrl/Shift-click) to overplot their data and fits. Clicking scans redraws
  from the last fit's outcomes and never refits.
* **Drop scans.** Remove Scan drops every selected scan. A scan removed from the
  project leaves the group too.

Ownership
---------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Who
     - Owns
   * - ``GroupFitWindow``
     - The group's scans, each scan's draft spec, the shared columns, and the
       last group fit's outcomes. This is display and draft state only; no
       model owns the group yet.
   * - ``GroupFitPresenter``
     - Nothing. It routes events to the window, and turns the Group Fit button
       into a ``GroupFitRequest`` built from the window's state.
   * - ``GroupFitModel``
     - Nothing persistent. It parses each spec with ``FitModel``'s translation,
       runs ``GroupFit`` and syncs the outcomes back.

Event flow
----------

.. mermaid::

    sequenceDiagram
        participant Tree as ProjectView
        participant PP as ProjectPresenter
        participant GP as GroupFitPresenter
        participant W as GroupFitWindow
        participant M as GroupFitModel (worker)

        Tree ->> PP: Add to Group Fit (uuid, name)*
        PP ->> GP: GroupScansEvent
        GP ->> W: add_scans_signal, then open
        Note over W: user edits the first scan, picks columns
        W ->> GP: Group Fit clicked
        GP ->> M: perform_group_fit(GroupFitRequest)
        M -->> GP: SyncGroupFitEvent(outcomes)
        GP ->> W: show_outcomes_signal (queued to the GUI thread)
        Note over W: selection changes redraw from cached outcomes

``RemoveRawScanEvent`` and ``RemoveProcessedScanEvent`` reach
``GroupFitPresenter`` too, which drops that scan from the window.

Suggest Params. follows the same shape:

.. mermaid::

    sequenceDiagram
        participant W as GroupFitWindow
        participant GP as GroupFitPresenter
        participant M as GroupFitModel (worker)

        W ->> GP: Suggest Params. clicked
        GP ->> M: suggest_peak_params(SuggestGroupPeakParamsRequest)
        Note over M: resolves the scan's columns, then FitModel.suggest_peak_params
        M -->> GP: SyncGroupPeakParamsEvent (or SyncGroupCosineParamsEvent)
        GP ->> W: peak_params_suggested_signal(uuid, ...)
        Note over W: applied only if that scan is still on screen

Key decisions
-------------

* **The first scan is the group model.** ``GroupFitModel`` turns a later scan's
  differences from the first scan into ``GroupFit`` overrides, so with chaining
  on everything else is seeded from the previous scan's fit. That only works if
  untouched fields really match the first scan, which is why first-scan edits
  flow to the other scans' unedited fields, and new scans are seeded from the
  first scan rather than from whichever scan is on screen.
* **A blank range is the whole scan.** Every scan in a group spans its own x,
  so no single starting range suits them all. The window starts with min/max
  empty, and ``GroupFitModel`` reads a blank bound as unbounded. The 1D Fitting
  tab is unchanged.
* **The current scan and the selection are separate.** The last scan clicked is
  the one the panel edits. Every selected scan is overplotted, each colored by
  its place in the list so it keeps its color as the selection changes.
* **Guesses have their own events.** ``FitModel`` publishes a guess as
  ``SyncPeakParamsEvent`` (and its cosine/background siblings), which the 1D
  Fitting tab fills its panel from. ``GroupFitModel`` reuses that code but
  republishes the result as ``SyncGroup*ParamsEvent``, subclasses that the
  broker dispatches separately because it matches on exact type. A guess asked
  for in one panel never lands in the other when both show the same scan.
* **The window caches outcomes, not the presenter.** They are display data for
  the canvas, like ``Plot1DView``'s fit curves. Re-selecting scans needs no
  event and no model call.
