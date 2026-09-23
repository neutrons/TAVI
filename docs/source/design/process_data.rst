ProcessOps
==========

.. class:: ProcessOps(tavi_data: TaviData, uuids: Sequence[UUID])

   Abstract base for operations that build derived :class:`ProcessedScan`
   objects from scans a project has already loaded. Lives in
   ``tavi.library.data.process_data``.

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
        ProcessOps <|-- AppendOp
        ProcessOps ..> TaviData : reads fetch_by_uuid()
        AppendOp ..> ProcessedScan : produces
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

``exec`` returns a ``ProcessedScan`` and registers it nowhere. Putting it in
``TaviData.processed_scans`` and announcing it is the model layer's job, because
the model layer owns ``TaviData`` — an operation only ever reads the pool it was
handed. A ``ProcessedScan`` may itself be an origin for the next combination,
which is why uuids resolve through ``fetch_by_uuid`` rather than against
``raw_scans`` alone.

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
