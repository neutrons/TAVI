Plot Data Model Philosophy
===========================

Overview
--------

A ``Plot`` is a *composition*, not a copy. It never owns measurement data —
it only records **which scans** to display and **which columns** of each
scan to use for x, y, and error bars. The actual numbers are pulled from the
referenced scan on demand, every time a render is needed.

This page documents the rule, why it exists, how multiple scans are overlaid
on one plot, and the intended path for how rebinning fits into this model
once that feature resumes.

Core Rule: Plot Holds No Data
------------------------------

``tavi.library.data.plot`` defines two classes:

.. code-block:: python

    class PlotSeries(BaseModel):
        """One scan's contribution to a Plot: which scan, and which columns of it to display."""

        source_scan_uuid: UUID
        scan_name: str
        normalized_by: Optional[str]
        normalized_by_value: Optional[float] = None
        x_name: str
        y_name: str
        error_name: str


    class Plot(BaseModel):
        """Composition of one or more PlotSeries displayed together. Holds no data itself."""

        uuid: UUID = UUIDFactory()
        series: list[PlotSeries]

Neither class carries an ``x``, ``y``, or ``err`` array. ``PlotSeries`` is a
pointer (``source_scan_uuid``) plus a column specification (``x_name`` /
``y_name`` / ``error_name`` / ``normalized_by``). ``Plot`` is nothing more
than a named list of these pointers.

.. mermaid::

    classDiagram
        class Plot {
            +UUID uuid
            +list~PlotSeries~ series
        }
        class PlotSeries {
            +UUID source_scan_uuid
            +str scan_name
            +Optional~str~ normalized_by
            +Optional~float~ normalized_by_value
            +str x_name
            +str y_name
            +str error_name
        }
        class Scan {
            <<RawScan today; ComboScan/ProcessedScan later>>
            +UUID uuid
            +ScanData data
            +TaviMetadata tavimeta
        }
        Plot "1" --> "1..*" PlotSeries : series
        PlotSeries "1" --> "1" Scan : source_scan_uuid

Why this matters:

- **No stale snapshots.** If a ``Plot`` cached arrays, those arrays would
  drift out of sync the moment the underlying scan or its selected columns
  changed. Pointing at the scan and re-reading columns at render time means
  there is only ever one source of truth.
- **No duplicated memory.** Scan data can be large; a plot that overlays the
  same scan on several axes should not carry N copies of it.
- **Multi-scan overlays fall out for free.** Because ``Plot.series`` is a
  list, adding a second scan to a plot is just appending a second
  ``PlotSeries`` — no special-casing required (see below).
- **Rebin/processing becomes a scan-level concern, not a plot-level one**
  (see `Future Direction`_).

Resolving a Series to Arrays
-----------------------------

Something still has to turn a ``PlotSeries`` into actual numbers before
matplotlib can draw them, but that must not become an excuse to hand the
presenter or view a live handle into a model's storage. Presenters are pure
orchestration and views are pure rendering — neither may own or reach into
state. So resolution is split across the event boundary instead of being a
method call across it:

``tavi.backend.model.plot_resolver`` holds two plain, stateless functions:

.. code-block:: python

    def resolve_series(series: PlotSeries, scans: dict[UUID, Scan]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Look up the scan this series points at and pull the x/y/err arrays named by it."""
        scan = scans[series.source_scan_uuid]
        x_name = "_" + series.x_name if series.x_name[0].isdigit() else series.x_name
        x = np.array(scan.data.data[x_name])
        y = np.array(scan.data.data[series.y_name])
        err = np.sqrt(np.abs(y))
        if series.normalized_by is not None:
            channel_data = np.array(scan.data.data[series.normalized_by])
            weight = series.normalized_by_value if series.normalized_by_value is not None else 1.0
            # Error propagation for the ratio y/channel_data, computed against the
            # pre-normalization y before it is overwritten below.
            err = weight * (err / channel_data)
            y = y * weight / channel_data
        return x, y, err


    def scans_for_plots(plots: list[Plot], scans: dict[UUID, Scan]) -> dict[UUID, Scan]:
        """Return the minimal {uuid: Scan} slice of `scans` referenced by any series in `plots`."""
        return {series.source_scan_uuid: scans[series.source_scan_uuid]
                for plot in plots for series in plot.series}

The ``x_name`` fixup mirrors the loader: column names beginning with a digit are
stored with a leading underscore so they remain attribute-accessible
(``2theta`` becomes ``_2theta``), and the series stores the user-facing name.

Whichever model actually owns scan storage — ``PlotModel`` (raw-scan-focus,
``update_fields``) or ``TaviProjectModel`` (re-focusing a previously saved
``Plot``) — calls ``scans_for_plots`` right before publishing, and attaches
the result to ``PlotFocusEvent.scans``:

.. code-block:: python

    self._event_broker.publish(PlotFocusEvent(plots=[plot], scans=scans_for_plots([plot], self._raw_scans)))

``EventBroker.publish`` deep-copies every event before delivering it (this
is not special-cased for ``scans`` — every event field already gets this),
so what the presenter receives is an independent snapshot, not a live
reference into ``self._raw_scans``.

``PlotterPresenter.handle_plot_focus`` then calls ``resolve_series`` once
per ``PlotSeries`` across ``e.plots``, reading only ``e.scans`` — never a
model, never a stored handle:

.. code-block:: python

    def handle_plot_focus(self, e: PlotFocusEvent) -> None:
        resolved = [(*resolve_series(series, e.scans), series) for plot in e.plots for series in plot.series]
        self._view.render_plots_signal.emit(resolved)

The view never sees a scan, a live model, or a ``PlotSeries`` from storage —
only the resolved arrays and the ``PlotSeries`` snapshot carrying the
display strings (``scan_name``, ``x_name``, ``y_name``, ``normalized_by``,
``error_name``) for that one render pass.

.. note::

   It is tempting to instead give the presenter a direct ``raw_scans``
   handle (mirroring how ``PlotModel`` is constructed with one), or to let
   the presenter call a model method synchronously for the resolved arrays.
   Both were tried and rejected: a handle makes the presenter a co-owner of
   model state, and a synchronous call-back defeats the async
   proxy/worker/event design (see ``tavi.meta.multithreading.proxy.Proxy``)
   that every other model method relies on. The event itself must already
   carry everything the presenter/view need.

Example: A Single Raw Scan
----------------------------

Focusing one ``RawScan`` produces a ``Plot`` with exactly one ``PlotSeries``:

.. code-block:: python

    scan = RawScan(
        uuid=UUID(value="scan-001"),
        data=ScanData(data={"qh": [1.0, 2.0, 3.0], "detector": [10.0, 40.0, 90.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=("qh", "detector"),
            normalization=("monitor", 1.0),
            friendly_name="scan 1",
            friendly_path="/exp1/scan0001.dat",
        ),
        prov=Provenance(raw_file="scan0001.dat", contributing_scans={UUID(value="scan-001"): 1}),
    )

    series = PlotSeries(
        source_scan_uuid=scan.uuid,
        scan_name=scan.tavimeta.friendly_name,
        normalized_by=scan.tavimeta.normalization[0],
        x_name="qh",
        y_name="detector",
        error_name="error",
    )
    plot = Plot(series=[series])

    x, y, err = resolve_series(plot.series[0], {scan.uuid: scan})
    assert list(x) == [1.0, 2.0, 3.0]

Example: Multiple Scans on One Plot
--------------------------------------

Overlaying scans is just adding another ``PlotSeries`` that points at a
different ``source_scan_uuid``. Nothing about ``Plot`` changes:

.. code-block:: python

    plot = Plot(
        series=[
            PlotSeries(
                source_scan_uuid=scan_1.uuid,
                scan_name=scan_1.tavimeta.friendly_name,
                normalized_by="monitor",
                x_name="qh",
                y_name="detector",
                error_name="error",
            ),
            PlotSeries(
                source_scan_uuid=scan_2.uuid,
                scan_name=scan_2.tavimeta.friendly_name,
                normalized_by="monitor",
                x_name="qh",
                y_name="detector",
                error_name="error",
            ),
        ]
    )

    assert len(plot.series) == 2
    assert plot.series[0].source_scan_uuid != plot.series[1].source_scan_uuid

    # The publishing model attaches both referenced scans to the event...
    scans = scans_for_plots([plot], {scan_1.uuid: scan_1, scan_2.uuid: scan_2})

    # ...and PlotterPresenter resolves + appends each series independently
    # against that snapshot, so both end up as separate error-bar containers
    # on the same axes.
    for series in plot.series:
        x, y, err = resolve_series(series, scans)
        view.append_plot(x, y, err, series.scan_name, series.normalized_by,
                          series.x_name, series.y_name, series.error_name)

Future Direction
------------------

The rebin controls in ``Plot1DView`` are currently disabled — rebinning
mutated ``x``/``y``/``err`` directly inside ``PlotModel.update_fields``
before ``Plot`` stopped owning data, and that path was removed rather than
carried forward half-adapted. Deciding what "rebin" means in a world of
scan-referencing plots is still open, but the shape of the answer is already
implied by the model:

**A rebin operation must not reach into ``Plot`` or ``PlotSeries`` — it
must produce a new scan for a ``PlotSeries`` to point at.**

Concretely: rebinning ``x``/``y``/``err`` from a source scan should create
a **``ProcessedScan``** — a scan-shaped object carrying the *rebinned*
columns, keyed by its own fresh ``UUID`` and recording the scan(s) it was
derived from (the existing ``Provenance.contributing_scans: dict[UUID, int]``
field already anticipates exactly this — a derived scan mapping back to one
or more contributing source uuids). ``ProcessedScan`` would sit in the same
role ``RawScan`` does today: immutable once produced, stored in a
scan-lookup dict keyed by uuid, and resolvable by ``resolve_series`` without
any change to that function — ``resolve_series``/``PlotFocusEvent.scans``
are already typed as ``Scan``, not ``RawScan``, precisely so a
``ProcessedScan`` slots in without touching either.

.. mermaid::

    classDiagram
        class PlotSeries {
            +UUID source_scan_uuid
        }
        class RawScan {
            +UUID uuid
            +ScanData data
        }
        class ProcessedScan {
            +UUID uuid
            +ScanData data
            +Provenance prov
        }
        PlotSeries "1" --> "1" RawScan : points at (today)
        PlotSeries "1" --> "1" ProcessedScan : points at (after rebin)
        ProcessedScan "1" --> "1..*" RawScan : derived from (prov.contributing_scans)

Sketch of the intended flow, once implemented:

.. code-block:: python

    # Not yet implemented — illustrates the intended shape of the feature.
    processed = rebin_tolerance(raw_scan, tolerance=0.02)

    assert isinstance(processed, ProcessedScan)
    assert processed.uuid != raw_scan.uuid
    assert raw_scan.uuid in processed.prov.contributing_scans

    # The PlotSeries is repointed at the processed scan — its column names
    # stay the same, only the referenced scan changes.
    rebinned_series = series.model_copy(update={"source_scan_uuid": processed.uuid})
    assert rebinned_series.source_scan_uuid == processed.uuid

    x, y, err = resolve_series(rebinned_series, {processed.uuid: processed})  # reads `processed`, not `raw_scan`

Under this design, "undo rebin" is just repointing the ``PlotSeries`` back
at the original ``source_scan_uuid`` — no data needs to be restored, because
the raw scan was never mutated in the first place.

Key Design Decisions
-----------------------

Plot is a view, scans are the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Plot``/``PlotSeries`` belong to the presentation/composition layer;
``RawScan`` (and, eventually, ``ProcessedScan``) belong to the data layer.
Turning one into the other is a pure function (``resolve_series``), never a
stateful lookup a presenter or view performs on its own.

Presenters and views hold no state — the event does
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only the models that own scan storage (``PlotModel``, ``TaviProjectModel``)
ever read ``raw_scans``. Whichever one publishes ``PlotFocusEvent`` attaches
exactly the scans its plots' series reference (``scans_for_plots``) to
``PlotFocusEvent.scans``. ``PlotterPresenter`` and ``Plot1DView`` resolve
against that per-event snapshot and nothing else — no live handle, no
synchronous call back into a model. This is the same reason
``RawScanFocusEvent`` already carries full ``RawScan`` objects rather than
uuids for the presenter to look up: the event is the only channel through
which rendering-layer code is allowed to see data.

Series list, not a single pointer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supporting overlays only required ``Plot.series`` to be a list instead of a
single field. No new plot subclass or "multi-plot" concept was needed.

Rebin output is a scan, not a plot mutation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Keeping ``Plot`` data-free only holds up if every transformation of the
underlying numbers (rebin included) is modeled as "produce a new scan,
point at it" rather than "give the plot new arrays." This is why rebin was
paused rather than re-wired directly onto ``PlotSeries`` when ``Plot`` lost
its arrays.
