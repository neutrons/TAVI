Selection → Visualization Flow
==============================

Overview
--------

Selecting an item in the project tree triggers a chain of typed events that
ultimately renders a plot in the ``Plot1DView`` widget.  The chain is fully
decoupled: no layer holds a direct reference to any other layer; each component
publishes an event and relies on the ``EventBroker`` singleton to deliver it
to the next stage.

Two variants exist depending on what was selected:

- **RawScan selected** — scan data must be converted into a plottable form
  before rendering (full chain)
- **Plot selected** — a pre-built plot object is retrieved directly and
  forwarded to the renderer, skipping the conversion step


Participants
------------

Five components take part, each with a single responsibility:

- **Project tree** — translates user interaction into a generic selection event
- **Load presenter** — collects the selected identifiers and publishes a focus event
- **Project model** — resolves each identifier to a typed domain object and
  re-publishes a type-specific event
- **Plot model** — builds a ``Plot`` (composition of ``PlotSeries``) from a
  focused raw scan (RawScan path only), and attaches the scans that
  ``Plot``'s series reference to the outgoing event
- **Plotter presenter + view** — resolves each series against the event's own
  scan snapshot, then clears the canvas and renders each resolved series

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

        User ->> TreeViewWidget: selects scan node (click or arrow keys)
        TreeViewWidget ->> LoadRawScanPresenter: selection signal (Qt)
        LoadRawScanPresenter ->> LoadRawScanPresenter: collect selected identifiers
        LoadRawScanPresenter ->> EventBroker: publish focus event

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> TaviProjectModel: resolve identifiers → RawScan objects
        TaviProjectModel ->> EventBroker: publish raw-scan focus event

        EventBroker ->> PlotModel: handle raw-scan focus event
        PlotModel ->> PlotModel: build Plot (PlotSeries pointing at scan.uuid)
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


Plot Selection — Short-Circuit Chain
-------------------------------------

When the selected identifier resolves to a pre-built plot object,
``TaviProjectModel`` publishes the plot focus event directly.
``PlotModel`` is never invoked — but ``TaviProjectModel`` still must gather
the scans referenced by that ``Plot``'s series (it owns ``raw_scans`` too)
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

        User ->> TreeViewWidget: selects plot node (click or arrow keys)
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
        PlotterPresenter ->> Plot1DView: append plot
        Plot1DView ->> Plot1DView: redraw canvas


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
keeps rendering concerns out of both the domain model and the view layer.

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

Mixed selection limitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user selects both scan and plot identifiers simultaneously, the project
model publishes both a raw-scan focus event and a plot focus event.  Because
the plot model reacts to the raw-scan event by publishing another plot focus
event, the plotter presenter receives two plot focus events in sequence and
clears the canvas between them.  Only the final group survives — a known
consequence of the clear-before-append design.
