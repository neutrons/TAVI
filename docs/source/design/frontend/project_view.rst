Project View ↔ NewRawScan Integration
====================================

Overview
--------

The Project View is updated reactively when new raw scans are loaded into the
project model. This interaction is mediated through an event-driven architecture
centered on the ``RawScanAppendEvent``.

The flow is:

1. Backend loads raw scans
2. Model emits ``RawScanAppendEvent``
3. Presenter listens and transforms the event
4. View updates the Project Tree

This decouples data ingestion from UI concerns.

Architecture
------------

The interaction spans three layers:

- **Model**: Owns raw scan data and emits events
- **Presenter**: Subscribes to events and updates the view
- **View (Project Tree)**: Displays scans hierarchically

Data Flow
---------

Loading raw scans triggers the following sequence:

.. code-block:: text

    FileMenuPresenter
        ↓
    TaviProjectModel.load_raw_scan_from_folder
        ↓
    RawScanLoadController.load_folder
        ↓
    (returns List[RawScan])
        ↓
    TaviProjectModel
        - stores scans in TaviData
        - emits RawScanAppendEvent per scan
        ↓
    EventBroker
        ↓
    LoadRawScanPresenter.update_treeview_data
        ↓
    ProjectView.add_raw_scan
        ↓
    TreeViewWidget.add_raw_scan

Model Responsibilities
----------------------

The ``TaviProjectModel`` is responsible for:

- Persisting scans in ``TaviData.raw_scans``
- Emitting a ``RawScanAppendEvent`` per scan

Each event contains:

- ``uuid``: Unique identifier
- ``friendly_name``: Display name
- ``friendly_path``: Logical grouping path

Implementation detail:

- Events are emitted *after* scans are stored
- One event is emitted per scan (not batched)

Event Definition
----------------

``RawScanAppendEvent`` is a typed model event:

.. code-block:: python

    class RawScanAppendEvent(Event):
        uuid: UUID
        friendly_name: str
        friendly_path: str

This event is the sole contract between backend data loading and UI updates.

Presenter Responsibilities
--------------------------

The ``LoadRawScanPresenter``:

- Registers with the ``EventBroker`` for ``RawScanAppendEvent``
- Maintains a local ``inventory`` of loaded scans
- Translates events into view updates

On event:

.. code-block:: python

    self._view.add_raw_scan(event.uuid,
                            event.friendly_name,
                            event.friendly_path)

The presenter does **no transformation** beyond forwarding metadata.

View Responsibilities (Project Tree)
------------------------------------

The Project View is implemented via ``TreeViewWidget`` and is responsible for:

- Maintaining hierarchical structure via ``path_map``
- Tracking nodes via ``uuid_map``
- Rendering scan entries under a ``/Raw`` root

Adding a Scan
~~~~~~~~~~~~~

All raw scans are inserted under a fixed root:

.. code-block:: text

    /Raw/<friendly_path>/<friendly_name>

The method:

.. code-block:: python

    add_raw_scan(uuid, name, path)

normalizes input by:

- Stripping leading ``/`` from ``path``
- Prefixing with ``/Raw``

Tree Construction
~~~~~~~~~~~~~~~~~

Tree nodes are created incrementally:

- Each path segment becomes a node if not already present
- Leaf nodes represent scans
- Duplicate names are allowed (UUID enforces uniqueness)

Key constraints:

- Duplicate UUID insertion raises an error
- Path structure is cached in ``path_map`` for O(1) lookup

Example
-------

Given:

- ``friendly_name = "scan1"``
- ``friendly_path = "/experiment/runA"``

The resulting tree:

.. code-block:: text

    /Raw
      └── experiment
            └── runA
                  └── scan1

Key Design Decisions
--------------------

Event-driven UI updates
~~~~~~~~~~~~~~~~~~~~~~~

The UI does not poll or directly query the model. Instead:

- The model *pushes* updates via events
- The presenter reacts asynchronously

This ensures:

- Loose coupling
- Easier extensibility (e.g., progress events, deletion events)

Per-scan event emission
~~~~~~~~~~~~~~~~~~~~~~~

Each scan generates its own event rather than batching.

Pros:

- Simpler UI updates
- Incremental rendering possible

Cons:

- Higher event volume for large loads

Strict UUID ownership
~~~~~~~~~~~~~~~~~~~~~

The view enforces uniqueness via ``uuid_map``:

- Prevents duplicate insertion bugs
- Enables future direct lookup (e.g., selection sync)

Summary
-------

The Project View does not directly depend on the loading mechanism. Instead, it
subscribes to ``RawScanAppendEvent`` and incrementally builds a hierarchical
representation of scans under the ``/Raw`` namespace.

This pattern cleanly separates:

- **Data ingestion** (controller + model)
- **Event propagation** (event broker)
- **Presentation logic** (presenter)
- **Rendering** (view)

----

Context Menu — Multi-Select Delete
===================================

Overview
--------

Right-clicking a scan node shows a context menu with a **Remove Item** action.
When multiple items are selected, the action removes all of them. When nothing
is selected, only the right-clicked item is removed.

Behaviour
---------

.. code-block:: text

    User selects one or more scan nodes (Ctrl+click / Shift+click)
        ↓
    User right-clicks any node
        ↓
    show_context_menu collects selected leaf indexes (column 0, has parent)
        ↓
    Falls back to right-clicked index if selection is empty
        ↓
    User clicks "Remove Item"
        ↓
    Each index converted to QPersistentModelIndex
        ↓
    remove_entry called per persistent index (invalid ones skipped)
        ↓
    uuid_map cleaned, tree rows removed

Selection Mode
--------------

``TreeViewWidget`` uses ``ExtendedSelection``, so users can multi-select with
``Ctrl+click`` or ``Shift+click`` before right-clicking.

Implementation Detail — Persistent Indexes
-------------------------------------------

Qt ``QModelIndex`` values are invalidated when rows are removed from the model.
Deleting item A shifts the row numbers of subsequent items, making previously
captured indexes stale.

To avoid this, every index is wrapped in a ``QPersistentModelIndex`` before the
deletion loop:

.. code-block:: python

    persistent = [QPersistentModelIndex(i) for i in to_delete]
    for pi in persistent:
        if pi.isValid():
            self.remove_entry(QModelIndex(pi))

``QPersistentModelIndex`` is kept valid by Qt as the model mutates. Any index
whose item was already deleted (e.g., a child whose parent was removed first)
returns ``isValid() == False`` and is silently skipped.

Recursive Delete
----------------

``remove_entry`` walks children depth-first before removing itself:

.. code-block:: python

    def remove_entry(self, index: QModelIndex) -> None:
        for row in range(self.treeModel.rowCount(index)):
            child = self.treeModel.index(row, 0, index)
            self.remove_entry(child)
        self._remove_index(index)

``_remove_index`` removes the row from the model and purges the UUID from
``uuid_map``.

Selecting a Parent and Child
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user selects a parent folder and a child scan, the parent's
``remove_entry`` call deletes all its children first. The child's persistent
index then becomes invalid and is skipped — no double-free or error.

Key Design Decisions
--------------------

Fallback to single-item delete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Right-clicking an item that is not part of the current selection deletes only
that item. This matches standard file-manager UX: a right-click that targets
an unselected item acts on that item alone, leaving the existing selection
unchanged.

Visibility of the menu action
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The "Remove Item" action appears only when there is at least one deletable
item (a leaf node that has a parent row). Right-clicking an empty area or a
root folder shows an empty menu.
