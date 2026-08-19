Data File View
===============

Overview
--------

``DataFileView`` is the panel that shows whichever scan is currently focused
by the user: its column data, a checklist of columns, and its metadata. The
view itself is agnostic to scan type — ``populate_columns``,
``populate_variables``, and ``populate_metadata`` take plain dicts/lists, not
a ``Scan`` object. It is driven entirely by ``DataFilePresenter``, which
subscribes to both ``RawScanFocusEvent`` (a scan focused directly from the
tree) and ``ActivePlotChangedEvent`` (a plot made active via the plotter's
"Current Plot" dropdown — see :doc:`visualization_flow`), and forwards
whichever scan is relevant to the view. The view holds no scan data of its
own between events.

Architecture
------------

- **Presenter** (``DataFilePresenter``): Subscribes to ``RawScanFocusEvent``
  and ``ActivePlotChangedEvent`` via the ``EventBroker`` and translates each
  event into view calls, through a shared ``_populate_from_scan`` helper.
  Holds no scan data itself.
- **View** (``DataFileView``): A ``QWidget`` with three regions — a data
  table, a variable checklist, and a tabbed metadata panel — plus a
  ``title_changed`` signal used to retitle its owning tab.

Data Flow
---------

.. code-block:: text

    RawScanFocusEvent(scans=[...])
        ↓
    EventBroker
        ↓
    DataFilePresenter.handle_raw_scan_focus
        - scans is empty  → view.clear_data(); view.set_title("Data File")
        - scans[0] is used, rest ignored
        ↓
    DataFilePresenter._populate_from_scan(scan)

    ActivePlotChangedEvent(plot, scans)
        ↓
    EventBroker
        ↓
    DataFilePresenter.handle_active_plot_changed
        - plot is None or has no series      → view.clear_data(); view.set_title("Data File")
        - plot.series[0].source_scan_uuid not in scans → view.clear_data(); view.set_title("Data File")
        - otherwise: look up that series' scan in `scans`
        ↓
    DataFilePresenter._populate_from_scan(scan)

    DataFilePresenter._populate_from_scan(scan):
        DataFileView.populate_columns(scan.data.data)
        DataFileView.populate_variables(list(scan.data.data.keys()))
        DataFileView.populate_metadata(scan.metadata.by_category())
        DataFileView.set_title(f"Data File ({scan.tavimeta.friendly_name})")
            ↓
        DataFileView.title_changed(title) --Qt signal--> TaviView relabels tab 0

Only one scan is ever displayed at a time — whichever scan the *active*
selection points at, regardless of whether that selection came from the
tree directly (``RawScanFocusEvent``) or from switching plots in the
plotter's dropdown while multiple scans/plots are focused
(``ActivePlotChangedEvent``).

Data Table and Variable Checklist
----------------------------------

``populate_columns`` rebuilds the data table with one column per key in
``scan.data.data``. ``populate_variables`` rebuilds the checklist with one
checked row per column name. Every cell in both tables has ``ItemIsEditable``
cleared — the panel is read-only.

Two interactions link the checklist back to the data table:

**Visibility.** Unchecking a row hides the matching data table column, matched by
header text, via ``_on_variable_check_changed``; rechecking shows it again. A
checkbox embedded in the checklist's corner button flips every row at once.

**Ordering.** Dragging a variable row to a new position moves the corresponding
data table column to the same position (``_on_variable_row_moved``). The two are
matched by *logical* index rather than visual position: rows and columns are
always populated together from the same ordered name list, so a row's logical
index is always the matching column's logical index regardless of how either
header has since been reordered.

Both are purely view concerns — neither affects the underlying data.

The **Restore Metadata From File** and **Save Modified Metadata To File** buttons
at the bottom of the panel exist but are disabled; metadata editing is not
implemented.

Populating the Metadata Widget
-------------------------------

The metadata panel renders one tab per **category**, and one table row per
field within that category. It is populated via:

.. code-block:: python

    DataFileView.populate_metadata(metadata: dict[str, dict])

``metadata`` must be ``{category_display_name: {field_name: value, ...}, ...}``.
Each top-level value must be a flat dict — anything else raises
``ValueError(f"Metadata tab {tab_name!r} must be a dict, not {type(fields).__name__!r}.")``.
Non-dict field values (e.g. a list) are rendered via ``str()`` in a single
table cell rather than expanded further.

When a focus event carries no scans, ``clear_data()`` empties both tables and
resets the metadata widget to a single tab labelled ``"Empty"``, which is also
the state the panel starts in.

Loaders and other code should not build this dict by hand. Instead, populate
``ScanMetadata`` correctly and call ``ScanMetadata.by_category()`` to get the
view-ready shape:

.. code-block:: python

    class ScanMetadata(BaseModel):
        data: Dict[str, Any]              # flat {field_name: value}
        categories: Dict[str, List[str]]  # display name -> list of `data` keys in that category

    def by_category(self) -> Dict[str, Any]:
        return {
            display_name: {key: self.data[key] for key in keys if key in self.data}
            for display_name, keys in self.categories.items()
        }

``data`` stays a flat ``{field_name: value}`` mapping — it is never nested by
category. ``categories`` only records, for each display name, which ``data``
keys should be grouped together for display. This is what decouples display
grouping from storage layout: a loader keeps a single flat dict of fields and
separately declares how those fields should be grouped into tabs.

Example — ``ORNLSpiceLoader.parse_metadata`` groups every field it parsed
under a single ``"ORNL Metadata"`` category:

.. code-block:: python

    data = metadata | {"errors": error_messages} | {"others": others}
    categories = {"ORNL Metadata": data.keys()}
    return ScanMetadata(data=data, categories=categories)

``scan.metadata.by_category()`` then returns
``{"ORNL Metadata": {...same flat dict, keyed by name...}}``, which
``DataFileView.populate_metadata`` renders as a single tab.

Flat Attribute Access Still Works
----------------------------------

Because ``data`` is never nested, ``ScanMetadata.__getattr__`` needs no
knowledge of ``categories`` at all — ``.``-delimited access (e.g.
``scan.metadata.def_x``) is a plain top-level lookup in ``data``. Callers
that just need a single field (loader internals, plotting code) never need to
know which category, if any, a field is grouped under for display.

Key Design Decisions
---------------------

Categories decouple storage from display
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Earlier, the metadata tab widget was built directly from the top-level keys
of ``ScanMetadata.data``, which forced storage layout and display grouping to
match exactly. ``categories`` + ``by_category()`` break that coupling: a
loader's internal data organization can change without changing how it's
displayed, and vice versa.

One scan at a time
~~~~~~~~~~~~~~~~~~~

``handle_raw_scan_focus`` only ever renders ``scans[0]``, and
``handle_active_plot_changed`` only ever renders the active plot's first
series' scan. The panel is a single-scan inspector, not a multi-scan
comparison view — even when the plotter is showing several plots at once,
the Data File tab always reflects exactly one of them (the active one).

Tab retitling is push-based, not polled
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DataFileView`` doesn't know it lives inside a ``QTabWidget`` — it just
emits ``title_changed`` on every repopulation (including the empty/cleared
case, which resets to plain ``"Data File"``). ``TaviView`` is the one that
knows the tab index and wires the signal to ``tabs.setTabText(0, title)``.
This keeps the view ignorant of its container, matching how it already
receives content purely through method calls rather than reaching out for
context.
