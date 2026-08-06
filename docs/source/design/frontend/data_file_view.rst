Data File View
===============

Overview
--------

``DataFileView`` is the panel that shows whichever scan is currently focused
by the user: its column data, a checklist of columns, and its metadata. The
view itself is agnostic to scan type — ``populate_columns``,
``populate_variables``, and ``populate_metadata`` take plain dicts/lists, not
a ``Scan`` object. It is driven entirely by ``DataFilePresenter``, which
today subscribes to ``RawScanFocusEvent`` (the only focus event currently
wired up) and forwards the focused scan's contents to the view. The view
holds no scan data of its own between events.

Architecture
------------

- **Presenter** (``DataFilePresenter``): Subscribes to ``RawScanFocusEvent``
  via the ``EventBroker`` and translates each event into view calls. Holds no
  scan data itself.
- **View** (``DataFileView``): A ``QWidget`` with three regions — a data
  table, a variable checklist, and a tabbed metadata panel.

Data Flow
---------

.. code-block:: text

    RawScanFocusEvent(scans=[...])
        ↓
    EventBroker
        ↓
    DataFilePresenter.handle_raw_scan_focus
        - scans is empty  → view.clear_data()
        - scans[0] is used, rest ignored
        ↓
    DataFileView.populate_columns(scan.data.data)
    DataFileView.populate_variables(list(scan.data.data.keys()))
    DataFileView.populate_metadata(scan.metadata.by_category())

Only the first scan in the event is displayed — the panel shows one scan at a
time, regardless of scan type.

Data Table and Variable Checklist
----------------------------------

``populate_columns`` rebuilds the data table with one column per key in
``scan.data.data``. ``populate_variables`` rebuilds the checklist with one
checked row per column name. Unchecking a row hides the matching data table
column (matched by header text) via ``_on_variable_check_changed``; rechecking
shows it again. Column visibility is purely a view concern — it does not
affect the underlying data.

Populating the Metadata Widget
-------------------------------

The metadata panel renders one tab per **category**, and one table row per
field within that category. It is populated via:

.. code-block:: python

    DataFileView.populate_metadata(metadata: dict[str, dict])

``metadata`` must be ``{category_display_name: {field_name: value, ...}, ...}``.
Each top-level value must be a flat dict — anything else raises ``ValueError``.
Non-dict field values (e.g. a list) are rendered via ``str()`` in a single
table cell rather than expanded further.

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

``handle_raw_scan_focus`` only ever renders ``scans[0]``. The panel is a
single-scan inspector, not a multi-scan comparison view.
