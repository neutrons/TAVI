Filter Model and Filtering the Project Browser
==============================================

Overview
--------

A filter narrows which raw scans the Project Browser lists. It is a named set of
criteria, each comparing one metadata entry of a scan (``key``) against a value
with an operator (``is``, ``is not``, ``has``, ``has not``, ``=``, ``!=``,
``>``, ``<``, ``>=``, ``<=``). A scan passes a filter when it satisfies every
criterion. Plots and fits are never filtered.

Saved filters live in ``TaviData.filters`` beside the scans, plots and fits, so
they are saved with the project like everything else, even though the project
tree never lists them. ``FilterModel`` owns them: it holds a live handle on
``TaviData.filters`` and saves and removes filters in it directly, and it tracks
which filter is selected and which is applied.

Data model
----------

.. mermaid::

    classDiagram
        class FilterOperator {
            <<StrEnum>>
            IS, IS_NOT, HAS, HAS_NOT
            EQUAL, NOT_EQUAL, GREATER, LESS, GREATER_EQUAL, LESS_EQUAL
            +is_numeric bool
        }
        class FilterCriterion {
            +key str
            +operator FilterOperator
            +value str
            +validate_criterion() list~str~
            +matches(metadata) bool
        }
        class ScanFilter {
            +uuid UUID
            +name str
            +criteria list~FilterCriterion~
            +validate_filter() list~str~
            +matches(scan) bool
        }
        class FilterModel {
            -_filters dict~UUID, ScanFilter~
            -_selected_uuid UUID
            -_applied_uuid UUID or None
            +save_filter(scan_filter)
            +remove_filter(filter_uuid)
            +select_filter(filter_uuid)
            +validate_filter(scan_filter)
            +apply_filter(scan_filter)
            +clear_filter()
        }
        class TaviData {
            +filters dict~UUID, ScanFilter~
            +filter_raw_scans(filters) list~UUID~
        }
        ScanFilter "1" *-- "*" FilterCriterion
        FilterCriterion --> FilterOperator
        TaviData "1" *-- "*" ScanFilter
        FilterModel ..> TaviData : edits filters

* Text operators (``is``, ``has`` and their negations) ignore case. Numeric
  operators compare both sides as numbers; a metadata value that is not a number
  does not match.
* A scan without the criterion's ``key`` never matches, whatever the operator,
  so ``is not`` does not sweep in every scan from an instrument that doesn't
  record the field.
* Every project starts with the read-only **Sample** filter (``scan >= 0``),
  selected and not applied. It is built into ``FilterModel`` and never stored in
  ``TaviData.filters``, so it cannot be saved over or removed, and is not written
  into project files.
* Sample's uuid, name and criterion are set once in ``application.yml`` under
  ``filter.sample``. ``sample_filter()`` builds the library's copy from it, and
  ``FilterView`` reads the same keys to list and lock its Sample entry. The view
  never imports the library's filter module, and the two can't drift apart.

Saving and removing
-------------------

.. mermaid::

    sequenceDiagram
        participant FP as FilterPresenter
        participant FM as FilterModel
        participant TD as TaviData.filters
        participant TPM as TaviProjectModel

        FP->>FM: save_filter(scan_filter)
        FM->>TD: filters[uuid] = scan_filter
        opt filter is new
            FM-->>FP: AddFilterEvent(filter)
        end

        FP->>FM: remove_filter(filter_uuid)
        opt it was applied
            FM->>TPM: ClearFilterEvent
        end
        Note over FM: falls back to selecting Sample if it was selected
        FM->>TD: del filters[uuid]
        FM-->>FP: RemoveFilterEvent(uuid)

Applying
--------

The filter panel asks ``FilterModel`` to act through its ``Proxy``: only that
model answers, so these are direct calls, not events. ``FilterModel`` then
broadcasts what is applied, because what the project lists is a concept other
models may come to depend on, the same way staging is broadcast.

.. mermaid::

    sequenceDiagram
        participant FP as FilterPresenter
        participant FM as FilterModel
        participant TPM as TaviProjectModel
        participant TD as TaviData
        participant LRS as ProjectPresenter
        participant PV as ProjectView

        FP->>FM: apply_filter(filter as shown)
        Note over FM: saves and selects it
        alt selected filter is invalid
            FM-->>FP: SyncFilterValidationEvent(errors)
        else valid
            opt a filter is already applied
                FM->>TPM: ClearFilterEvent
                Note over TPM,PV: relists every raw scan, as below
            end
            FM->>TPM: ApplyFilterEvent(filter)
            TPM->>TPM: ClearFocusEvent
            TPM->>LRS: ClearRawScansEvent
            LRS->>PV: clear_raw_scans()
            TPM->>TD: filter_raw_scans([filter])
            loop each passing raw scan
                TPM->>LRS: AddRawScanEvent(uuid, name, path)
                LRS->>PV: add_raw_scan(uuid, name, path)
            end
        end

Key decisions
-------------

**Apply only adds.** ``ApplyFilterEvent`` follows the lexicon's additive rule:
``TaviProjectModel`` keeps every filter applied since the last
``ClearFilterEvent`` and lists the scans that pass all of them. ``FilterModel``
applies one filter at a time, so it publishes a clear before replacing the
applied filter, the same way a new selection is ``ClearFocusEvent`` followed by
a focus.

**The model owns what is listed; the tree only mirrors it.** ``TaviProjectModel``
decides which raw scans are listed, asking ``TaviData.filter_raw_scans`` which
pass. The tree holds no filter state of its own: on a filter change it is told to
drop every raw scan (``ClearRawScansEvent``) and is then sent an
``AddRawScanEvent`` for each one that passes, exactly as if they had just been
loaded. The scans that fail stay in ``TaviData``; nothing announces
``RemoveRawScanEvent``, which would tell ``PlotModel`` and the fit presenters
that the scan had left the project.

**Filtering clears the focus.** Someone who filters wants to stop looking at what
the filter excludes, and the tree's selection goes with its rows, so every
relisting starts with ``ClearFocusEvent``. The tree drops its selection without
publishing a ``FocusEvent`` of its own, so the chain does not re-enter itself.

**Loading while filtered.** ``TaviProjectModel`` remembers the applied filters
it was told about. A load stores every new scan but announces only those that
pass; the rest are listed once the filters that exclude them are cleared.

**Filters live in the project, and FilterModel edits them.** There is no reason
to split them out of the main datastore just because the tree doesn't show them:
keeping them in ``TaviData`` saves them with the project. No other model reads or
writes them, so ``FilterModel`` edits its handle on ``TaviData.filters`` directly
rather than routing saves and removals through ``TaviProjectModel``, and
announces them itself with ``AddFilterEvent`` and ``RemoveFilterEvent``. They
stay out of ``TaviData.purge`` and ``remove_items``, which serve the tree.

**Validation is reported, not raised.** Proxy calls run on a worker thread and
return nothing to the caller, so ``validate_filter`` publishes its outcome as
``SyncFilterValidationEvent``. ``apply_filter`` validates first and publishes the
same event, applying nothing, when the selected filter is invalid.

The filter panel
----------------

``FilterView`` holds no filters of its own. Its New, Del, Apply and Clear buttons
and its dropdown only report what the user asked for; ``FilterPresenter`` turns
each into a ``FilterModel`` call, and the panel changes only when the model's
events come back:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - User action
     - ``FilterPresenter`` calls
     - Panel updates on
   * - New
     - ``save_filter`` with a blank criterion
     - ``AddFilterEvent``: listed and picked
   * - Pick in dropdown
     - ``select_filter``
     - (nothing to show)
   * - Apply
     - ``apply_filter`` with the filter as the panel shows it
     - ``ApplyFilterEvent``: named as applied; or ``SyncFilterValidationEvent``: errors shown
   * - Clear
     - ``clear_filter``
     - ``ClearFilterEvent``: "Default Applied", greyed out
   * - Del
     - ``remove_filter``
     - ``RemoveFilterEvent``: unlisted, Sample picked

These events are published from the model's worker thread, so the presenter
emits the view's ``*_signal`` signals rather than calling it, and Qt runs the
update on the GUI thread. The same goes for the project tree: relisting runs on
``FilterModel``'s worker thread, so every ``ProjectView`` change (add, remove,
clear) is a signal, and the queued calls run in the order they were emitted.

**Apply carries the panel's edits.** Every Proxy call runs on its own worker
thread, so a ``save_filter`` followed by ``apply_filter`` could apply before the
save landed. ``apply_filter`` therefore takes the filter as the panel shows it,
and saves, selects, validates and applies it in one call. An invalid filter is
still saved and selected, so the user's edits survive while they fix it. For
Sample the model applies its own built-in copy, whatever the panel sends.
