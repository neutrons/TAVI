.. _tavifields:

=======================
Fields and Validation
=======================

This page documents the user-editable fields of the TAVI GUI, where their values
come from, and where each is validated.

Plotter Controls
================

The 1D Plotter panel (``Plot1DView``) is the only panel with editable fields
today. Its controls are collected into a
:class:`tavi.library.data.plot.PlotFields` snapshot by ``get_plot_fields()`` and
passed to ``PlotModel.update_fields``. There is no "plot" button — a field
losing focus emits ``fields_focus_changed``, which triggers the round trip.

.. list-table:: Axis Fields
  :header-rows: 1

  * - Field
    - Type
    - Value Origin
    - Default
    - Validation
  * - X Axis
    - String
    - column name of the focused scan
    - ``scan.tavimeta.default_axis[0]``
    - must be a key of ``scan.data.data``
  * - Y Axis
    - String
    - column name of the focused scan
    - ``scan.tavimeta.default_axis[1]``
    - must be a key of ``scan.data.data``

.. list-table:: Rebin Fields
  :header-rows: 1

  * - Field
    - Type
    - Value Origin
    - Default
    - Validation
  * - Rebin Mode
    - Enum
    - radio group -> :class:`tavi.library.data.enum.rebin_mode.RebinMode`
      (``none`` / ``tolerance`` / ``equal_step``)
    - ``No Rebin`` (``RebinMode.NONE``)
    - constrained by the radio group
  * - Rebin Start
    - String
    -
    - ``"0"``
    - not consumed yet
  * - Rebin Stop
    - String
    -
    - ``"2"``
    - not consumed yet
  * - Rebin Step
    - String
    -
    - ``"0.02"``
    - not consumed yet

.. note::

   The rebin controls are currently **disabled**. Rebinning previously mutated
   the plotted arrays in place; that path was removed when ``Plot`` stopped
   owning data. The intended replacement — producing a new scan for the series to
   point at — is described in :doc:`frontend/plot_data_model`.

.. list-table:: Preset (Normalization) Fields
  :header-rows: 1

  * - Field
    - Type
    - Value Origin
    - Default
    - Validation
  * - Preset Type
    - Enum
    - :class:`tavi.library.data.enum.preset_type.PresetType`
      (``none`` / ``normalize``)
    - ``none``
    - constrained by the combo box
  * - Preset Channel
    - String
    - dropdown repopulated from the focused scan's column names
    - first column of the focused scan
    - must be a key of ``scan.data.data`` when the type is ``normalize``
  * - Preset Value
    - String
    -
    - ``"1"``
    - must parse as ``float`` when the type is ``normalize``

Where Validation Happens
========================

TAVI does **not** use Qt validators on these fields. All of them are free-form
``QLineEdit`` / ``QComboBox`` values; validation is a backend concern, performed
in :meth:`tavi.backend.model.plot_model.PlotModel._resolve_series_update`:

* ``x_axis`` and ``y_axis`` are stripped and looked up in ``scan.data.data``.
* When ``preset_type`` is ``NORMALIZE``, ``preset_channel`` must also be a column
  and ``preset_value`` must parse as a float.

If any check fails the update is rejected wholesale: ``_apply_fields_to_plot``
returns ``None``, no ``PlotFocusEvent`` is published, and the canvas keeps showing
the last valid plot. An invalid entry is therefore a no-op rather than an error
dialog — the fields resync to the plotted values on the next successful render.

Field Synchronization
=====================

Whatever is actually rendered is written back into the controls, so the fields
never drift from the plot:

``sync_axis_fields(x_name, y_name)``
    Called per rendered series from ``Plot1DView._render_plots``. Reflects the
    plotted column names in the axis fields.

``sync_preset_fields(normalized_by, normalized_by_value)``
    Also called from ``_render_plots``. Sets the preset type to ``NORMALIZE``
    when the series carries a normalization channel and ``NONE`` otherwise, then
    fills in the channel and value.

``set_preset_channel_options(columns)``
    Called by ``PlotterPresenter.handle_raw_scan_focus``. Repopulates the channel
    dropdown from the newly focused scan's columns.

``reset_controls_to_defaults()``
    Also called on raw-scan focus, before the dropdown is repopulated. Returns the
    rebin radio group and the preset type/value to their defaults.

All four block widget signals while writing, so a programmatic sync never
re-triggers ``fields_focus_changed`` and loops back into the model.

Filter Panel
============

``FilterView`` renders Temperature, a ± tolerance, and an IPTS number, plus
**Clear Filters** and **Apply Filters** buttons. The widgets exist and hold
default values (3, 0.01, 1234), but nothing is wired to a presenter or model yet
— the panel is a layout placeholder.

Data File Panel
===============

``DataFileView`` has no editable fields. Its data table and metadata tables are
explicitly non-editable (``ItemIsEditable`` cleared on every item). Two
interactions do exist:

* **Variable checklist** — unchecking a row hides the matching data table column,
  matched by header text. Purely a view concern; it does not affect the data.
* **Variable row reordering** — dragging a row in the variable table moves the
  corresponding data table column to the same position.

The **Restore Metadata From File** and **Save Modified Metadata To File** buttons
are present but disabled; metadata editing is not implemented.
