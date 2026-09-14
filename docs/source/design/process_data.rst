ProcessData
===========

.. class:: ProcessData(tavi_data: TaviData)

   Builds derived :class:`ProcessedScan` objects from scans a project has
   already loaded. Lives in ``tavi.library.data.process_data``.

   :param tavi_data: The project data holding the raw and processed scans that
      uuids are resolved against, via ``TaviData.fetch_by_uuid``.
   :type tavi_data: TaviData

Overview
--------

``ProcessData.append(uuids, columns)`` appends several scans into one. It takes
the origin scans **by uuid**, in the order their rows should be appended, and the
**column names** to carry over; the result is a ``ProcessedScan`` whose every
requested column is the origins' columns of that name laid end to end.

.. code-block:: python

    processor = ProcessData(tavi_data)
    combined = processor.append([scan_a.uuid, scan_b.uuid], ["qh", "en", "detector"])

    assert combined.data.qh == list(scan_a.data.qh) + list(scan_b.data.qh)
    assert combined.prov.contributing_scans == {scan_a.uuid: 1, scan_b.uuid: 1}
    assert combined.tavimeta.friendly_name == f"{scan_a_name}+{scan_b_name}"

This is the "append" case of the scan-level processing anticipated in
:doc:`frontend/plot_data_model` — a derived scan a ``PlotSeries`` can point at,
carrying its own fresh ``UUID`` and recording where it came from, with no raw
scan mutated and nothing added to ``Plot``/``PlotSeries``.

.. mermaid::

    classDiagram
        class ProcessData {
            +TaviData tavi_data
            +append(uuids, columns) ProcessedScan
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
        ProcessData ..> TaviData : reads fetch_by_uuid()
        ProcessData ..> ProcessedScan : produces
        ProcessedScan "1" --> "1..*" Scan : derived from (prov.contributing_scans)

Public Methods
--------------

.. method:: append(uuids: Sequence[UUID], columns: Sequence[str]) -> ProcessedScan

   Append several scans into one, column by column.

   The origin scans are taken in the order given, and each requested column of
   the result is the origins' columns of that name laid end to end. Values are
   stored as plain floats, since loaders may hand back numpy arrays. Columns not
   named in ``columns`` are dropped.

   :param uuids: Origin scan uuids, in the order their rows should be appended.
   :type uuids: Sequence[UUID]
   :param columns: Names of the columns to carry over. At least two, the first
      two becoming the result's ``tavimeta.default_axis``.
   :type columns: Sequence[str]
   :returns: A ``ProcessedScan`` with a fresh uuid, each origin uuid recorded in
      ``prov.contributing_scans`` at weight 1, and the origins' friendly names
      joined with ``+`` as its ``tavimeta.friendly_name``. Its ``prov.raw_file``
      and ``tavimeta.friendly_path`` are empty.
   :rtype: ProcessedScan
   :raises ValueError: If fewer than two columns are requested.
   :raises KeyError: If a uuid is not present in the data pool.

Key Design Decisions
--------------------

Producing is not storing
~~~~~~~~~~~~~~~~~~~~~~~~

``append`` returns a ``ProcessedScan`` and registers it nowhere. Putting it in
``TaviData.processed_scans`` and announcing it is the model layer's job, because
the model layer owns ``TaviData`` — ``ProcessData`` only ever reads the pool it
was handed. A ``ProcessedScan`` may itself be an origin for the next
combination, which is why uuids resolve through ``fetch_by_uuid`` rather than
against ``raw_scans`` alone.

Columns are the caller's responsibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A requested column missing from an origin contributes nothing rather than
raising, and a column no origin carries comes back empty. This leaves the
result's columns at differing lengths, which misaligns rows that are meant to be
read together — but which columns belong together is a question only the caller
can answer, so ``append`` combines what it is asked for and does not second-guess
the selection. Listing a uuid twice is likewise permitted: its rows appear twice
while ``prov.contributing_scans`` keeps a single entry for it.

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
