Data File View
===============

Overview
--------

``DataFileView`` is the panel that shows the raw scan currently focused by the
user: its column data, a checklist of columns, and its metadata. It is driven
entirely by ``DataFilePresenter``, which subscribes to ``RawScanFocusEvent``
and forwards the focused scan's contents to the view. The view holds no scan
data of its own between events.

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
time.

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
        data: Dict[str, Any]        # however the loader wants to organize fields
        categories: Dict[str, str]  # display name -> key in `data` for that category

    def by_category(self) -> Dict[str, Any]:
        return {display_name: self.data.get(data_key, {}) for display_name, data_key in self.categories.items()}

``categories`` maps a front-end display name to the key in ``data`` holding
that category's fields. This is what decouples display grouping from storage
layout: a loader is free to organize ``data`` however it likes (flat, nested
by instrument section, etc.) and simply declares which top-level keys should
become tabs, and what they're called.

Example — ``ORNLSpiceLoader.parse_metadata`` groups everything under a single
``"ORNL Metadata"`` category:

.. code-block:: python

    data = {"ORNL Metadata": metadata}
    categories = {"ORNL Metadata": "ORNL Metadata"}
    return ScanMetadata(data=data, categories=categories)

``scan.metadata.by_category()`` then returns
``{"ORNL Metadata": {...the flat metadata dict...}}``, which
``DataFileView.populate_metadata`` renders as a single tab.

Flat Attribute Access Still Works
----------------------------------

``ScanMetadata.__getattr__`` supports ``.``-delimited access (e.g.
``scan.metadata.def_x``) regardless of how fields are grouped into
categories: it checks top-level ``data`` first, then searches inside each
category's dict. Callers that just need a single field (loader internals,
plotting code) never need to know which category a field lives under.

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
