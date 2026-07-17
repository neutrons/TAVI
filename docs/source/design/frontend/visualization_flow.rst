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
- **Plot model** — converts raw scan domain objects into renderable plot objects
  and re-publishes for the renderer (RawScan path only)
- **Plotter presenter + view** — clears the canvas and renders each plot object


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

        User ->> TreeViewWidget: clicks scan node
        TreeViewWidget ->> LoadRawScanPresenter: selection signal (Qt)
        LoadRawScanPresenter ->> LoadRawScanPresenter: collect selected identifiers
        LoadRawScanPresenter ->> EventBroker: publish focus event

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> TaviProjectModel: resolve identifiers → RawScan objects
        TaviProjectModel ->> EventBroker: publish raw-scan focus event

        EventBroker ->> PlotModel: handle raw-scan focus event
        PlotModel ->> PlotModel: extract axes data from scan
        PlotModel ->> PlotModel: build plot objects
        PlotModel ->> EventBroker: publish plot focus event

        EventBroker ->> PlotterPresenter: handle plot focus event
        PlotterPresenter ->> Plot1DView: clear canvas
        loop for each plot object
            PlotterPresenter ->> Plot1DView: append plot
        end
        Plot1DView ->> Plot1DView: redraw canvas


Plot Selection — Short-Circuit Chain
-------------------------------------

When the selected identifier resolves to a pre-built plot object,
``TaviProjectModel`` publishes the plot focus event directly.
``PlotModel`` is never invoked.

.. mermaid::

    sequenceDiagram
        participant User
        participant TreeViewWidget
        participant LoadRawScanPresenter
        participant EventBroker
        participant TaviProjectModel
        participant PlotterPresenter
        participant Plot1DView

        User ->> TreeViewWidget: clicks plot node
        TreeViewWidget ->> LoadRawScanPresenter: selection signal (Qt)
        LoadRawScanPresenter ->> EventBroker: publish focus event

        EventBroker ->> TaviProjectModel: handle focus event
        TaviProjectModel ->> TaviProjectModel: resolve identifiers → Plot objects
        TaviProjectModel ->> EventBroker: publish plot focus event

        EventBroker ->> PlotterPresenter: handle plot focus event
        PlotterPresenter ->> Plot1DView: clear canvas
        PlotterPresenter ->> Plot1DView: append plot
        Plot1DView ->> Plot1DView: redraw canvas


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

``PlotModel`` exists to convert raw scan domain objects into plot
view-model objects.  This keeps rendering concerns out of both the domain
model and the view layer.

Clear-before-append semantics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The plotter presenter always clears the canvas before appending.  Every plot
focus event represents a *replacement*, not an *addition*.  Overplot behaviour
(accumulating multiple series on the same axes) requires a separate code path.

Mixed selection limitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user selects both scan and plot identifiers simultaneously, the project
model publishes both a raw-scan focus event and a plot focus event.  Because
the plot model reacts to the raw-scan event by publishing another plot focus
event, the plotter presenter receives two plot focus events in sequence and
clears the canvas between them.  Only the final group survives — a known
consequence of the clear-before-append design.
