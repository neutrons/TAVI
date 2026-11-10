ProcessOps
==========

.. class:: ProcessOps(tavi_data: TaviData, uuids: Sequence[UUID])

   Abstract base for operations that build derived :class:`ProcessedScan`
   objects from scans a project has already loaded. Lives in
   ``tavi.library.data.process_op``; each concrete operation is a module of its
   own under ``tavi.library.data.algorithms`` — ``append_ops.AppendOp``,
   ``equal_step_bin.EqualStepRebinOps`` and ``bin_by_tol.ToleranceRebinOps``.
   ``tavi.library.data.combine.CombineOps`` sits over the three as the single
   entry point a caller reaches for, resolving which of them the given
   parameters describe and storing what comes back.

   :param tavi_data: The project data holding the raw and processed scans that
      uuids are resolved against, via ``TaviData.fetch_by_uuid``.
   :type tavi_data: TaviData
   :param uuids: Origin scan uuids, in the order the operation takes them.
   :type uuids: Sequence[UUID]

Overview
--------

An operation is constructed with everything it needs, then run: ``exec()``
returns a ``ProcessedScan`` and registers it nowhere. ``AppendOp`` is the first
such operation — it appends several scans into one, taking the origin scans **by
uuid**, in the order their rows should be appended, and the **column names** to
carry over; the result is a ``ProcessedScan`` whose every requested column is the
origins' columns of that name laid end to end.

.. code-block:: python

    combined = AppendOp(tavi_data, [scan_a.uuid, scan_b.uuid], ["qh", "en", "detector"]).exec()

    assert combined.data.qh == list(scan_a.data.qh) + list(scan_b.data.qh)
    assert combined.prov.contributing_scans == {scan_a.uuid: 1, scan_b.uuid: 1}
    assert combined.tavimeta.friendly_name == f"{scan_a_name}+{scan_b_name}"

This is the "append" case of the scan-level processing anticipated in
:doc:`frontend/plot_data_model` — a derived scan a ``PlotSeries`` can point at,
carrying its own fresh ``UUID`` and recording where it came from, with no raw
scan mutated and nothing added to ``Plot``/``PlotSeries``.

.. mermaid::

    classDiagram
        class ProcessOps {
            <<abstract>>
            +TaviData tavi_data
            +Sequence~UUID~ uuids
            +int size
            +validate()*
            +exec()* ProcessedScan
        }
        class AppendOp {
            +Sequence~str~ columns
            +validate()
            +exec() ProcessedScan
        }
        class TaviData {
            +dict~UUID,RawScan~ raw_scans
            +dict~UUID,ProcessedScan~ processed_scans
            +fetch_by_uuid(uuid) Scan
        }
        class ProcessedScan {
            +UUID uuid
            +ScanData data
            +ScanMetadata metadata
            +TaviMetadata tavimeta
            +Provenance prov
        }
        class CombineOps {
            +Sequence~str~ columns
            +RebinMode mode
            +validate()
            +exec() ProcessedScan
        }
        ProcessOps <|-- AppendOp
        ProcessOps <|-- CombineOps
        ProcessOps ..> TaviData : reads fetch_by_uuid()
        AppendOp ..> ProcessedScan : produces
        CombineOps ..> AppendOp : delegates to
        CombineOps ..> TaviData : writes processed_scans
        ProcessedScan "1" --> "1..*" Scan : derived from (prov.contributing_scans)

ProcessOps Members
------------------

.. attribute:: size

   Number of origin uuids the operation was given, for GUI generation.

.. method:: validate() -> None

   Raise if the operation's inputs cannot produce a scan. Subclasses define
   what "usable" means for their own inputs.

.. method:: exec() -> ProcessedScan

   Run the operation and return the derived scan, registered nowhere.

AppendOp
--------

.. class:: AppendOp(tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str])

   Append several scans into one, column by column.

   The origin scans are taken in the order given, and each requested column of
   the result is the origins' columns of that name laid end to end. Values are
   stored as plain floats, since loaders may hand back numpy arrays. Columns not
   named in ``columns`` are dropped.

   :param tavi_data: The pool the origin uuids are looked up in.
   :type tavi_data: TaviData
   :param uuids: Origin scan uuids, in the order their rows should be appended.
   :type uuids: Sequence[UUID]
   :param columns: Names of the columns to carry over. At least two, the first
      two becoming the result's ``tavimeta.default_axis``.
   :type columns: Sequence[str]

   .. method:: validate() -> None

      :raises ValueError: If fewer than two columns are requested.
      :raises KeyError: If a uuid is not present in the data pool, or an origin
         does not carry one of the requested columns.

   .. method:: exec() -> ProcessedScan

      Calls ``validate()`` first, then lays the origins' columns end to end.

      :returns: A ``ProcessedScan`` with a fresh uuid, each origin uuid recorded
         in ``prov.contributing_scans`` at weight 1, and the origins' friendly
         names joined with ``+`` as its ``tavimeta.friendly_name``. Its
         ``prov.raw_file`` and ``tavimeta.friendly_path`` are empty.
      :rtype: ProcessedScan
      :raises ValueError: If fewer than two columns are requested.
      :raises KeyError: If a uuid is not present in the data pool, or an origin
         does not carry one of the requested columns.

CombineOps
----------

.. class:: CombineOps(tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str], norm: str = "", tol: Optional[float] = None, step: Optional[float] = None, start: Optional[float] = None, finish: Optional[float] = None, multiplication: float = 1.0)

   Combine one or more raw scans into a single ``ProcessedScan``, registered in
   ``tavi_data.processed_scans``. Lives in ``tavi.library.data.combine``.

   Which combination is run is read off the parameters rather than named, so a
   caller describes the result it wants instead of picking a class:

   .. list-table::
      :header-rows: 1
      :widths: 20 25 55

      * - Given
        - ``mode``
        - Delegate
      * - ``tol``
        - ``RebinMode.TOLERANCE``
        - ``ToleranceRebinOps(tavi_data, uuids, columns, norm, tol, multiplication)``
      * - ``step``
        - ``RebinMode.EQUAL_STEP``
        - ``EqualStepRebinOps(tavi_data, uuids, columns, norm, step, start, finish, multiplication)``
      * - neither
        - ``RebinMode.NONE``
        - ``AppendOp(tavi_data, uuids, columns)``

   ``CombineOps`` holds no binning maths of its own. It resolves the mode,
   refuses parameters that mode has no use for, delegates, and stores.

   .. code-block:: python

       result = CombineOps(
           tavi_data,
           [scan_a.uuid, scan_b.uuid],
           ["qh", "detector", "monitor"],
           norm="monitor",
           step=0.1,
           multiplication=10000,
       ).exec()

       assert tavi_data.fetch_by_uuid(result.uuid) is result
       assert result.tavimeta.normalization == ("monitor", 10000)

   :param tavi_data: The pool the origin uuids are looked up in, **and** the pool
      the result is written into.
   :type tavi_data: TaviData
   :param uuids: Origin scan uuids, at least one, each of which must be a ``RawScan``.
   :type uuids: Sequence[UUID]
   :param columns: Names of the columns to carry over, at least two.
   :type columns: Sequence[str]
   :param norm: Column each bin is weighted by. Required by both rebins, and
      refused in append mode, which has no bin to weigh.
   :type norm: str
   :param tol: Selects the tolerance rebin; how close two ``columns[0]`` values
      must be to count as the same point. Mutually exclusive with ``step``.
   :type tol: Optional[float]
   :param step: Selects the equal-step rebin; the bin width.
   :type step: Optional[float]
   :param start: Outer edge of the first equal-step bin; read off the data if left out.
   :type start: Optional[float]
   :param finish: Outer edge of the last equal-step bin; read off the data if left out.
   :type finish: Optional[float]
   :param multiplication: Scale applied to the per-``norm`` y column and its
      ``error``. Refused in append mode, which has nothing to scale.
   :type multiplication: float

   .. attribute:: mode

      The ``RebinMode`` the parameters resolve to, so a GUI can ask a built
      operation what it will do before running it. Raises ``ValueError`` if both
      ``tol`` and ``step`` were given.

   .. attribute:: columns

      The columns the result carries, the first being the axis binned on.

   .. method:: validate() -> None

      Check that each parameter belongs to the mode it resolves to, then defer to
      the delegate for everything else — positive, finite ``tol``/``step``/
      ``multiplication``, ``finish > start``, ``norm`` among ``columns`` and not
      the y column, no carried ``error`` column, and every origin carrying every
      requested column. None of those checks is duplicated here.

      :raises ValueError: If ``tol`` and ``step`` are both given; if no uuids were
         given; if an origin is not a ``RawScan``; if ``start`` or ``finish`` was
         given without a ``step``; or if ``norm`` or a ``multiplication`` other
         than 1 was given with no rebin to apply it.
      :raises KeyError: If a uuid is not present in the data pool, or an origin
         does not carry one of the requested columns.

   .. method:: exec() -> ProcessedScan

      Calls ``validate()``, runs the delegate, and writes the result into
      ``tavi_data.processed_scans`` under its fresh uuid before returning it. A
      refused combination stores nothing.

      :returns: The delegate's ``ProcessedScan`` unaltered — its columns, its
         ``error`` column where a rebin computed one, its
         ``prov.contributing_scans`` and its friendly name.
      :rtype: ProcessedScan

CombineModel
------------

.. class:: CombineModel(tavi_data: TaviData)

   The backend model a GUI combination goes through. Lives in
   ``tavi.backend.model.combine_model``, behind ``CombineModelInterface`` and
   its ``CombineModelProxy``, and is constructed in ``tavi.__main__`` with
   ``TaviProjectModel.get_tavi_data_handle()`` — the whole project store rather
   than one handle into it, because a combination resolves its origins across
   every store and writes its result back into one.

   ``CombineModel`` holds no combining logic. Its single call,
   ``combine(request)``, translates either side of ``CombineOps``:

   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Step
        - What happens
      * - Columns
        - ``[x_axis, y_axis]``, then the normalization column if it names a
          third one. The first two become the result's ``default_axis``, and
          the append drops every column not listed — so a ``norm`` left out
          would be gone before a rebin could weigh by it.
      * - Numbers
        - The request carries them as the text the user typed.
          ``rebin_mode`` decides which are read: a ``tolerance``, or a
          ``step`` with optional ``start``/``finish``. A value left in the
          other mode's fields is never parsed, so it cannot reach
          ``CombineOps`` as a second mode.
      * - Failure
        - Every ``ValueError``/``KeyError`` the ops raise describes what the
          user asked for, so each is published as a ``ReportErrorEvent`` and
          returned as ``ModelResponse(ERROR, message)`` rather than raised
          into a worker thread nothing is watching.
      * - Success
        - ``CombineOps.exec`` has already registered the result in
          ``tavi_data.processed_scans``; the model publishes
          ``AddProcessedScanEvent(uuid, friendly_name, friendly_path)``.

   .. code-block:: python

       CombineModel(tavi_data).combine(
           CombineRequest(
               source_scan_uuids=[scan_a.uuid, scan_b.uuid],
               x_axis="qh",
               y_axis="detector",
               normalization="monitor",
               rebin_mode=RebinMode.EQUAL_STEP,
               step="0.02",
               start="0",
               finish="2",
           )
       )

   ``CombineRequest`` (``tavi.library.data.combine``) is the request itself —
   what the user asked for, as they expressed it. It subscribes to nothing:
   running a combination is a request exactly one model answers, so it is a
   direct call through the proxy, and the result is broadcast as an event. See
   :doc:`/guides/event_broker`.

From the Combine tab
~~~~~~~~~~~~~~~~~~~~

``CombinePresenter`` (``tavi.frontend.presenter.combine_presenter``) owns the
``CombineView`` and is the only thing holding the two vocabularies together:
the view works in scan **names**, the model in **uuids**.

.. mermaid::

    sequenceDiagram
        actor User
        participant Tree as ProjectPresenter
        participant Broker as EventBroker
        participant Panel as CombinePresenter
        participant Model as CombineModel
        participant Ops as CombineOps

        User->>Tree: select raw scans
        Tree->>Broker: publish FocusEvent(ids)
        Broker-->>Panel: FocusRawScanEvent(scans)
        Note over Panel: remembers uuid → name,<br/>offers the scans' columns,<br/>seeds the axis fields
        User->>Panel: Add Scans
        Note over Panel: stages the focused scans;<br/>the view lists their names
        User->>Panel: Combine
        Panel->>Model: combine(CombineRequest)
        Model->>Ops: exec()
        Ops-->>Model: ProcessedScan, registered in TaviData
        Model->>Broker: publish AddProcessedScanEvent
        Broker-->>Tree: AddProcessedScanEvent
        Tree->>Tree: add_combined() — a row under /Combined
        User->>Tree: select the Combined row
        Tree->>Broker: publish FocusEvent(ids)
        Broker-->>Broker: TaviProjectModel routes it as FocusProcessedScanEvent
        Note over Panel: not subscribed — a processed scan<br/>cannot be staged to combine again

Selecting a row under ``/Combined`` focuses the processed scan the same way a
raw one is focused: ``TaviProjectModel`` publishes ``FocusProcessedScanEvent``,
``PlotModel`` answers with a preview plot of its ``default_axis`` — the
combination's x and y columns — and from there it is staged, edited, saved and
shown in the data tab like any raw scan (see
:doc:`frontend/visualization_flow`). It travels on its own event rather than
``FocusRawScanEvent`` because ``CombinePresenter`` subscribes to that one, and
``CombineOps`` refuses processed origins.

The staged list is the user's own basket, so it outlives a selection: a
``ClearFocusEvent`` forgets only which scans *Add Scans* would stage next. A
``RemoveRawScanEvent`` does unstage its scan, the list never naming data
nothing can resolve. The axis fields are re-seeded from each focused scan's
``tavimeta.default_axis``, as the plotter's are, so the panel cannot be left
naming the previous experiment's columns.

**Append** sends ``RebinMode.NONE`` with no normalization and no scale. The
panel's Normalization and Multiplication fields describe a rebin's weighting,
and ``CombineOps`` refuses a weighting an append cannot honour rather than
dropping it — so the presenter does not offer it one. The default
normalization is ``mcu``.

Key Design Decisions
--------------------

One class per operation
~~~~~~~~~~~~~~~~~~~~~~~

Each operation is its own ``ProcessOps`` subclass rather than a method on a
shared processor class. The base holds only what every operation needs — the
pool and the origin uuids — while operation-specific inputs (``columns`` for
append, a bin width for rebin) are constructor arguments of the subclass. An
operation is therefore a value that can be built, inspected (``size``) and
validated before it is run, which is what a GUI needs in order to offer the
operation before the user commits to it.

``validate`` is not enforced by the base
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``exec`` calls ``validate`` itself; the base class does not wrap it in a
template method. Operations differ in when their inputs can be checked — some
cheaply up front, some only once the origins are resolved — so forcing a single
validate-then-run order into the base would constrain subclasses without
protecting callers, who are free to call ``validate()`` on its own beforehand.

Producing is not storing
~~~~~~~~~~~~~~~~~~~~~~~~

``exec`` returns a ``ProcessedScan`` and registers it nowhere. An operation only
ever reads the pool it was handed, because the model layer owns ``TaviData``. A
``ProcessedScan`` may itself be an origin for the next combination, which is why
uuids resolve through ``fetch_by_uuid`` rather than against ``raw_scans`` alone.

``CombineOps`` is the one exception, and deliberately so. It is not another
primitive but the entry point a caller reaches for when it wants a combined scan
to *exist* in the project, and the rule above is what had left
``TaviData.processed_scans`` written by nobody: every op returned a scan and
every caller dropped it. Pushing the write into the one operation whose whole
purpose is to produce a stored result keeps that store from depending on each
caller remembering, while the primitives it delegates to stay pure and
separately testable. Announcing the new scan remains the model layer's job —
``CombineOps`` writes the dict, it publishes no event.

A separate exception, on a different axis: the rebinning operations
``EqualStepRebinOps`` and ``ToleranceRebinOps`` reduce each bin by summing raw
counts against the monitor they were counted with, which is only meaningful on
data nothing has already reduced, so they refuse any origin that is not a
``RawScan``. ``CombineOps`` applies that check in every mode, including the
append one, so that one façade means raw scans whichever operation it resolves
to. Chaining is still available in the other direction — a rebinned scan is an
ordinary ``ProcessedScan`` and may be appended directly through ``AppendOp``,
plotted or fitted like any other.

The mode is inferred, and inapplicable parameters are refused
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CombineOps`` takes no mode argument: ``tol`` means a tolerance rebin, ``step``
an equal-step one, and neither a plain append. The parameters already say which
combination is meant, and a separate mode argument would only create a second
place for them to disagree — a ``RebinMode.TOLERANCE`` carrying a ``step`` has no
reading. ``mode`` is still reported as a ``RebinMode`` so a caller can inspect a
built operation.

Having inferred it, ``CombineOps`` raises on any parameter that mode cannot use
rather than ignoring it: ``start`` or ``finish`` without a ``step``, or ``norm``
or a ``multiplication`` other than 1 with no rebin to apply them. Silently
dropping a ``norm`` would hand back an unnormalized scan to a caller who asked
for a normalized one, with nothing in the result saying so.

Every origin must carry every requested column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``validate()`` resolves each origin and raises ``KeyError`` — naming the scan,
the missing columns and the columns it does carry — if any requested column is
absent. Skipping the missing column instead would leave the result's columns at
differing lengths, silently misaligning rows that are meant to be read together;
a scan whose ``qh`` is two values longer than its ``detector`` is not something a
caller can detect after the fact, so the operation refuses up front rather than
producing it. Because the check needs the origins resolved, it lives in
``validate()`` alongside the column-count check rather than in ``__init__``, and
a GUI can call ``validate()`` on its own to offer the operation before the user
commits to it.

Listing a uuid twice is still permitted: its rows appear twice while
``prov.contributing_scans`` keeps a single entry for it.

At least two columns, for the default axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The origins' ``default_axis`` need not survive the combination, since the axis
columns may not have been requested. The first two requested columns are used
instead — they are guaranteed to be present in the result — which is why fewer
than two columns raises rather than leaving a derived scan the frontend cannot
plot.

Only what describes the combination is kept
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ScanMetadata`` is not carried over, and ``tavimeta.normalization`` is left
unset. These fields (scan number, temperature, monitor channel, ...) describe an
individual origin rather than the combination; carrying one origin's value
forward would be a lie. The origins remain recoverable per-uuid from
``TaviData`` through ``prov.contributing_scans``.

Provenance is the origin scans
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``prov.contributing_scans`` maps each origin uuid to weight 1 — the field's
intended use for derived scans, appending never modifying weights. ``raw_file``
is empty, because a combined scan has no file of its own. The origins' friendly
names appear in the combined scan's ``tavimeta.friendly_name``, joined with
``+``.
