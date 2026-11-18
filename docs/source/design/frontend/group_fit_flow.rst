Group Fit Window
================

Overview
--------

A **group fit** fits a sweep of scans (typically in temperature or field) in one
go, with :ref:`GroupFit <group_fit>` doing the fitting. It replaces the
deprecated sequential fit described in :doc:`fit_data_model`: it has its own
window, its own model and presenter, and the fit itself is display only. The one thing it saves into
``TaviData`` is a *summary*: one fitted parameter against a column's per-scan
average, listed under Fits.

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
  from the last fit's outcomes and never refits. The panel shows the current
  scan's own fitted values, uncertainties and reduced chi-squared, as the 1D
  Fitting tab does after a fit.
  Hide Label, beside the plot's toolbar, drops the legend, as in the 1D
  plotter.
* **Drop scans.** Remove Scan drops every selected scan. A scan removed from the
  project leaves the group too.
* **Summarize.** The line just above the Group Fit button reads "Summarize fit
  results: <parameter> vs <column> [Generate]". The parameter drop-down offers
  what the last group fit produced (``peak1 center``, ``peak1 FWHM``,
  ``peak1 Int. Inten.``, background terms, reduced chi-squared), and the column
  is any scan column, e.g. ``temp``. Generate saves a ``FitSummary`` under Fits
  with one point per fitted scan: the parameter with its 1-sigma, against the
  column's mean over the whole scan with its standard deviation. Selecting it in
  the tree draws it in the 1D plotter, y vs x with error bars on both axes;
  several selected summaries are overlaid, one color each. The Data File tab
  lists the first one's points as columns: scan, x, x err, y, y err.

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
       runs ``GroupFit`` and syncs the outcomes back. For a summary it averages
       the column from each scan and asks for the result to be saved.
   * - ``TaviProjectModel``
     - Saved summaries, as ``FitSummary`` entries in ``TaviData.fits``.

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

Generate turns the window's cached outcomes into a summary saved in the project:

.. mermaid::

    sequenceDiagram
        participant W as GroupFitWindow
        participant GP as GroupFitPresenter
        participant M as GroupFitModel (worker)
        participant TP as TaviProjectModel
        participant PP as ProjectPresenter

        W ->> GP: Generate clicked
        GP ->> M: summarize_group_fit(SummarizeGroupFitRequest)
        Note over M: each scan's parameter, and the column's mean/std from its data
        M ->> TP: SaveFitSummaryEvent(FitSummary)
        TP ->> PP: AddFitEvent
        Note over PP: listed under /Fits

Selecting a summary in the tree draws it in the 1D plotter and lists it in the
Data File tab:

.. mermaid::

    sequenceDiagram
        participant PP as ProjectPresenter
        participant TP as TaviProjectModel
        participant PL as PlotterPresenter
        participant V as Plot1DView
        participant DF as DataFilePresenter

        PP ->> TP: FocusEvent(ids)
        TP ->> PL: ClearFocusEvent
        TP ->> PL: FocusFitSummaryEvent(summaries)
        PL ->> V: add_fit_summaries_signal
        Note over V: y vs x with x and y error bars, alone on the canvas
        TP ->> DF: FocusFitSummaryEvent(summaries)
        Note over DF: first summary's FitSummary.columns() as the data table

``RemoveFitEvent`` for a shown summary takes it off the canvas and empties the
Data File tab again.

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
* **A summary holds its own values.** A ``FitEntry`` points at its scans and is
  recomputed when selected. A ``FitSummary`` is a table of numbers, like a
  combined scan holds its own data: nothing to recompute, and it outlives the
  scans it came from. It shares ``TaviData.fits`` with fit entries so it is
  listed under Fits, and code that reads members checks for ``FitEntry`` first.
* **The request carries the results.** Generate sends the last group fit's
  ``FitResultSummary`` per scan, taken from the window's outcomes, the same way
  Group Fit sends the window's specs. ``GroupFitModel`` stays stateless and only
  reads the averaged column from the scans. The drop-down's labels are the keys
  of ``FitResultSummary.parameters()``, so view and model name a parameter the
  same way.
* **A summary takes the canvas alone.** Its axes are a fitted parameter and a
  column's average, which share no scale with a scan's series, so like a 2D map
  it is drawn in place of any 1D series focused with it. It has its own
  ``FocusFitSummaryEvent`` rather than riding on ``FocusFitEvent``: there is no
  curve to recompute and no series for ``PlotModel`` to stage.
* **Fitted values on display are not starting values.** The panel shows the
  current scan's fit result in the same value fields that hold its starting
  values. Leaving a scan stores only the fields the user changed since it was
  shown (``_store_current`` compares against what was loaded). If the fitted
  values were stored as the spec, they would differ from the first scan's, so
  they would turn into per-scan overrides and switch chaining off for those
  parameters. Group Fit therefore always reruns from the starting values the
  user set.
