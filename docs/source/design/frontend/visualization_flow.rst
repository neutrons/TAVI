Selection → Visualization Flow
==============================

Overview
--------

Selecting an item in the project tree triggers a chain of typed events that
ultimately renders a plot in the ``Plot1DView`` widget.  The chain is fully
decoupled: no layer holds a direct reference to any other layer; each component
publishes an event and relies on the ``EventBroker`` singleton to deliver it
to the next stage.

The tree supports multi-select: every currently-selected scan or plot is
focused together, in the order the user selected them (see
`Selection order is preserved across multi-select`_ below).

Two variants exist depending on what was selected:

- **RawScan selected** — each selected scan's data must be converted into a
  plottable form before rendering (full chain)
- **Plot selected** — pre-built plot objects are retrieved directly and
  forwarded to the renderer, skipping the conversion step

Whichever batch is focused, the plotter's "Current Plot" dropdown lists every
plot in it, and exactly one is "active" at a time — see
`Switching the Active Plot`_ for how that selection is changed without
re-resolving or re-rendering the whole batch.


Participants
------------

Six components take part, each with a single responsibility:

- **Project tree** — translates user interaction into a generic selection event,
  preserving the order items were selected in
- **Load presenter** — collects the selected identifiers and publishes a focus event
- **Project model** — resolves each identifier to a typed domain object and
  re-publishes a type-specific event
- **Plot model** — builds one unsaved preview ``Plot`` (composition of
  ``PlotSeries``, keyed by that scan's own uuid) per focused raw scan
  (RawScan path only), and attaches the scans that those ``Plot``\(s)'
  series reference to the outgoing event
- **Plotter presenter + view** — resolves each series against the event's own
  scan snapshot, clears the canvas and renders each resolved series, and
  tracks which of the focused plots is "active"
- **Data file presenter + view** — mirrors whichever plot is currently
  "active" into the Data File tab, retitled to that plot's source scan

Neither the project model nor the plot model ever hands a live
``raw_scans``/``plots`` handle to a presenter. Whatever a presenter or view
needs to render travels *inside* the event that triggers it — see
:doc:`plot_data_model` for why.


RawScan Selection — Full Chain
-------------------------------

.. mermaid::

    sequenceDiagram
        participant User
        participant TreeViewWidget
        participant LoadRawScanPresenter
        participant EventBroker
        participant TaviProjectModel
        participant PlotModel
        participant PlotterPresenter
        participant Plot1DView

        User ->> TreeViewWidget: selects one or more scan nodes (click, ctrl/shift-click, arrow keys)
        TreeViewWidget ->> LoadRawScanPresenter: selection signal (Qt)
        LoadRawScanPresenter ->> LoadRawScanPresenter: collect selected identifiers, in selection order
        LoadRawScanPresenter ->> EventBroker: publish focus event

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> TaviProjectModel: resolve identifiers → RawScan objects
        TaviProjectModel ->> EventBroker: publish raw-scan focus event

        EventBroker ->> PlotModel: handle raw-scan focus event
        loop for each focused scan
            PlotModel ->> PlotModel: build one single-series preview Plot (PlotSeries pointing at scan.uuid)
        end
        PlotModel ->> PlotModel: gather referenced scans (scans_for_plots)
        PlotModel ->> EventBroker: publish plot focus event (plots + scans)

        EventBroker ->> PlotterPresenter: handle plot focus event
        loop for each series across plots
            PlotterPresenter ->> PlotterPresenter: resolve_series(series, event.scans)
        end
        PlotterPresenter ->> Plot1DView: clear canvas
        loop for each resolved series
            PlotterPresenter ->> Plot1DView: append plot
        end
        Plot1DView ->> Plot1DView: redraw canvas
        PlotterPresenter ->> Plot1DView: set_plot_options(labels, default_index)
        PlotterPresenter ->> EventBroker: publish ActivePlotChangedEvent(scan)
        EventBroker ->> DataFilePresenter: handle active plot changed


Plot Selection — Short-Circuit Chain
-------------------------------------

When the selected identifier(s) resolve to pre-built plot objects,
``TaviProjectModel`` publishes the plot focus event directly.
``PlotModel`` is never invoked — but ``TaviProjectModel`` still must gather
the scans referenced by those ``Plot``\(s)' series (it owns ``raw_scans`` too)
before publishing, for the same reason ``PlotModel`` does on the other path.

.. mermaid::

    sequenceDiagram
        participant User
        participant TreeViewWidget
        participant LoadRawScanPresenter
        participant EventBroker
        participant TaviProjectModel
        participant PlotterPresenter
        participant Plot1DView

        User ->> TreeViewWidget: selects one or more plot nodes (click, ctrl/shift-click, arrow keys)
        TreeViewWidget ->> LoadRawScanPresenter: selection signal (Qt)
        LoadRawScanPresenter ->> EventBroker: publish focus event

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> TaviProjectModel: resolve identifiers → Plot objects
        TaviProjectModel ->> TaviProjectModel: gather referenced scans (scans_for_plots)
        TaviProjectModel ->> EventBroker: publish plot focus event (plots + scans)

        EventBroker ->> PlotterPresenter: handle plot focus event
        loop for each series across plots
            PlotterPresenter ->> PlotterPresenter: resolve_series(series, event.scans)
        end
        PlotterPresenter ->> Plot1DView: clear canvas
        loop for each resolved series
            PlotterPresenter ->> Plot1DView: append plot
        end
        Plot1DView ->> Plot1DView: redraw canvas
        PlotterPresenter ->> Plot1DView: set_plot_options(labels, default_index)
        PlotterPresenter ->> EventBroker: publish ActivePlotChangedEvent(scan)


Switching the Active Plot
--------------------------

Both chains above end by rendering every focused plot and picking one as
"active" — the one whose source scan drives the Data File tab. Picking a
*different* already-focused plot from the "Current Plot" dropdown does not
replay either chain: nothing needs re-resolving, since every plot in the
dropdown was already resolved when the batch was focused. Only the one
newly-active plot needs to be looked up and announced.

The presenter publishes a single ``FocusActivePlotEvent(uuid)`` — it does not
need to know (and ``PlotFocusEvent`` carries no flag saying) which model
produced the currently-focused batch:

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant EventBroker
        participant PlotModel
        participant TaviProjectModel
        participant DataFilePresenter

        User ->> Plot1DView: picks a different "Current Plot" entry
        Plot1DView ->> PlotterPresenter: plot_combo_index_changed(index)
        PlotterPresenter ->> PlotterPresenter: _active_plot_uuid = _focused_plot_uuids[index]
        PlotterPresenter ->> EventBroker: publish FocusActivePlotEvent(uuid)

        EventBroker ->> PlotModel: handle active-plot focus event
        EventBroker ->> TaviProjectModel: handle active-plot focus event
        Note over PlotModel,TaviProjectModel: whichever one owns uuid resolves it;<br/>the other no-ops
        Note over PlotModel,TaviProjectModel: first_contributing_scan(plot, scans) — plot.series[0]'s source scan
        Note over PlotModel,TaviProjectModel: publish ActivePlotChangedEvent(scan)

        EventBroker ->> DataFilePresenter: handle active plot changed
        DataFilePresenter ->> DataFilePresenter: populate from scan

``FocusActivePlotEvent`` is registered by *both* ``TaviProjectModel`` and
``PlotModel`` — each checks membership first (``e.uuid in self.tavi_data.plots``,
or a linear scan over ``self._last_plots``) and returns without publishing if
the uuid isn't one of its own, rather than indexing straight in and letting a
miss raise. This works because saved-plot uuids (freshly minted via
``uuid4()`` on save) and preview-plot uuids (borrowed from the ``RawScan``
they preview — see `Preview plots are keyed by their source scan's uuid`_
below) never collide: exactly one model ever recognizes a given uuid as its
own, so the presenter never has to guess which model to ask, and neither
model needs to reason about the other's storage to decline gracefully.

No ``PlotFocusEvent`` is published on this path — the canvas is left alone.
Only ``ActivePlotChangedEvent`` fires, which is why switching the active plot
never re-clears or re-renders ``Plot1DView``.


Add Plot — Saving a New Plot Entry
------------------------------------

Clicking **Add Plot** on the plotter panel runs the reverse direction of the
chains above: instead of a selection producing a rendered plot, the
*currently rendered* plot is captured and persisted as a new, independent
entry in the project.

.. mermaid::

    sequenceDiagram
        participant User
        participant Plot1DView
        participant PlotterPresenter
        participant EventBroker
        participant TaviProjectModel
        participant LoadRawScanPresenter
        participant ProjectView

        User ->> Plot1DView: clicks "Add Plot"
        Plot1DView ->> PlotterPresenter: plot_clicked signal
        PlotterPresenter ->> PlotterPresenter: deep-copy _current_plot under a fresh uuid
        PlotterPresenter ->> EventBroker: publish SavePlotEvent(plot)

        EventBroker ->> TaviProjectModel: handle save-plot event
        TaviProjectModel ->> TaviProjectModel: store plot in TaviData.plots
        TaviProjectModel ->> TaviProjectModel: friendly_name = "_".join(series.scan_name) + "_Plot"
        TaviProjectModel ->> EventBroker: publish PlotAppendEvent(uuid, friendly_name, friendly_path="")

        EventBroker ->> LoadRawScanPresenter: handle plot append event
        LoadRawScanPresenter ->> ProjectView: add_plot(uuid, friendly_name, friendly_path)

The new plot is then selectable like any other — its own future
``FocusEvent`` will resolve it through the "Plot Selection" chain above.

Why a fresh uuid on every click
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``_current_plot`` may already be a *saved* plot (the user re-selected an
existing ``Plots`` tree entry, then clicked **Add Plot** again).
``TaviProjectModel.tavi_data.plots`` is keyed by uuid, and ``ProjectView``
raises if a uuid is inserted twice — reusing the uuid would silently
overwrite the original entry (model) or crash the tree (view) on the second
click. Deep-copying under a new uuid before publishing means every click
produces an independent entry, regardless of whether the source plot was
ephemeral (fresh off a raw-scan selection) or already saved.

``friendly_path`` is always empty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Saved plots have no folder/grouping requirement yet, so
``TaviProjectModel`` always publishes ``friendly_path=""``. ``ProjectView``
still places every plot under the ``Plots`` tree root regardless.


Key Design Decisions
--------------------

Typed event routing in the project model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The project model dispatches to different downstream events based on the
runtime type of each resolved object.  This keeps the upstream trigger generic
— the caller does not need to know whether it selected a scan or a plot —
while allowing downstream consumers to specialize.

PlotModel as an adapter
~~~~~~~~~~~~~~~~~~~~~~~~

``PlotModel`` exists to convert raw scan domain objects into ``Plot``
compositions (``Plot``/``PlotSeries`` — see :doc:`plot_data_model`).  This
keeps rendering concerns out of both the domain model and the view layer. A
multi-scan focus produces one single-series preview ``Plot`` per scan (not
one multi-series ``Plot``), so each run stays independently selectable in
the "Current Plot" dropdown and independently resolvable by uuid.

Preview plots are keyed by their source scan's uuid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A preview ``Plot`` built from a raw scan is deliberately never written into
``TaviData.plots`` — only clicking **Add Plot** persists it (see below), and
until then it exists only in ``PlotModel._last_plots``. Rather than minting
a fresh, TaviData-unaware uuid for it (``PlotSeries`` and ``UUIDFactory``
would happily do this), ``_preview_plot_for_scan`` gives the preview
``Plot`` the *same* uuid as the one ``RawScan`` it previews — a preview is a
1:1, single-series wrapper around one run, so this loses no information and
makes the uuid trivially real: it is always the uuid of something already
sitting in ``TaviData.raw_scans``.

This is what lets ``FocusActivePlotEvent`` (see `Switching the Active Plot`_)
be handled by both models with neither reasoning about the other: since a
preview plot's uuid is a real ``TaviData.raw_scans`` key,
``TaviProjectModel._handle_active_plot_focus_event`` cannot use
``fetch_by_uuid`` here — it would resolve to the ``RawScan``, the wrong type,
instead of telling us the uuid isn't a plot at all. It checks
``e.uuid in self.tavi_data.plots`` directly instead, and no-ops on a miss.
``PlotModel`` mirrors this by scanning its own ``_last_plots`` for a matching
uuid. Neither model imports, mentions, or reasons about the other's storage;
each just knows its own identity space.

Events carry their own data; presenters/views hold none
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PlotFocusEvent`` carries both the focused ``Plot``\(s) and the ``Scan``
objects their series reference (``scans``), gathered by whichever model
publishes it. ``PlotterPresenter`` resolves each series against that
snapshot and forwards resolved arrays to the view. Neither the presenter nor
the view ever holds a live handle to ``raw_scans``, nor calls a model
synchronously to fetch data — every value they need arrives already inside
the event that triggered them.

Clear-before-append semantics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The plotter presenter always clears the canvas before appending.  Every plot
focus event represents a *replacement*, not an *addition*.  Overplot behaviour
(accumulating multiple series on the same axes) requires a separate code path.

Controls are reset on scan focus, synced on plot render
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two focus paths touch the plotter controls differently, and the ordering
matters because a RawScan focus is always immediately followed by a plot focus
(``PlotModel`` reacts to the first by publishing the second).

- **RawScan focus** — ``PlotterPresenter.handle_raw_scan_focus`` calls
  ``reset_controls_to_defaults()`` (rebin radio back to *No Rebin*, rebin
  start/stop/step back to ``0``/``2``/``0.02``, preset type back to ``NONE``,
  preset value back to ``"1"``), then ``set_preset_channel_options(...)`` to
  repopulate the channel dropdown from the newly focused scan's columns. It does
  **not** derive normalization from ``scan.tavimeta.normalization``.
- **Plot focus** — ``Plot1DView._render_plots`` calls ``sync_axis_fields(x_name,
  y_name)`` and ``sync_preset_fields(normalized_by, normalized_by_value)`` per
  rendered series, so the axis and preset controls end up reflecting what is
  actually on the canvas.

So a scan carrying a normalization channel does briefly show preset type
``NONE`` — the reset — before the plot focus event lands and
``sync_preset_fields`` promotes it to ``NORMALIZE``. Because both events are
dispatched synchronously through the broker, the user only ever sees the final
state.

Both sync methods and the reset block widget signals while writing, so
programmatically updating a control never re-emits ``fields_focus_changed`` and
loops back into the model.

.. note::

   ``_render_plots`` calls both sync methods **once per series**, so with a
   multi-series overlay the controls end up showing the last series' values.
   The controls describe a single series; the overlay case has no defined
   answer yet.

Selection order is preserved across multi-select
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``TreeViewWidget`` tracks a ``_selection_order`` list alongside Qt's own
selection model, appending a uuid when it's selected and removing it when
it's deselected. ``get_selected_items()`` returns items in that order rather
than Qt's (unordered) ``selectedIndexes()`` order, and drops any uuid no
longer selected (e.g. a removed tree entry) before returning. Deselecting an
item and reselecting it moves it to the end — the tracked order reflects the
current *click* order, not first-ever-selected order. This determines the
order plots appear in in the "Current Plot" dropdown and, on a genuinely new
selection, which one defaults to active (index 0 — see
``handle_plot_focus``'s ``new_uuids[0]`` fallback).

The data widget follows the active plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DataFilePresenter`` subscribes to ``ActivePlotChangedEvent`` in addition to
``RawScanFocusEvent``. In a multi-select scenario this is what keeps the Data
File tab in sync with whichever plot the user has made active via the
dropdown, rather than freezing on whatever scan was focused first.
``DataFileView.set_title`` emits a ``title_changed`` signal that
``TaviView`` connects to retitle the tab itself (e.g. ``"Data File
(scan_name)"``); a ``None`` scan (nothing active) resets the title to plain
``"Data File"``. See :doc:`data_file_view` for the presenter/view details.

``ActivePlotChangedEvent`` carries a scan, not a plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The data widget only ever displays data that lives on a ``Scan`` — it has no
use for a ``Plot`` object itself. So rather than carry the active ``Plot``
plus a ``scans`` snapshot for ``DataFilePresenter`` to resolve against (which
also meant handling the case where that ``Plot`` is an unsaved preview with
nowhere persistent to live), every publisher of ``ActivePlotChangedEvent``
resolves ``first_contributing_scan(plot, scans)`` (``plot_resolver.py``) —
the scan backing the plot's first series — itself, and the event carries
that ``Scan`` directly (``scan: Optional[Scan]``). A ``Plot`` always has at
least one series pointing at a real scan, whether the plot is a saved
``TaviData`` entry or an unsaved preview, so this always resolves; ``None``
only ever means no plot is currently active. ``DataFilePresenter`` is left
with a single branch: ``None`` clears the view, otherwise populate from the
scan.

Mixed selection limitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user selects both scan and plot identifiers simultaneously, the project
model publishes both a raw-scan focus event and a plot focus event.  Because
the plot model reacts to the raw-scan event by publishing another plot focus
event, the plotter presenter receives two plot focus events in sequence and
clears the canvas between them.  Only the final group survives — a known
consequence of the clear-before-append design. Selection-order tracking does
not change this: it orders items *within* one focus event, not across the
two separately-dispatched events here.
