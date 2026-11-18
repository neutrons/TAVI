.. _group_fit:

+++++++++
Group Fit
+++++++++

``GroupFit`` processes a group of scans (typically a sweep in temperature or field) in one call. For
every scan it fits several y columns against one shared x column and reduces other columns to a
per-scan average. Each single fit is delegated to :class:`tavi.library.fit.fit.Fit` (see
:ref:`fit_package`), so models are described exactly as ``Fit.fit`` takes them.

Source module:

.. code-block:: text

   tavi.library.fit.group_fit


Minimal Example
===============

.. code-block:: python

   from tavi.library.fit import GroupFit, GroupFitParams, ModelName

   params = GroupFitParams(
       x_col="s1",  # None -> each scan's tavimeta.default_axis[0]
       fits={
           "detector": [
               (ModelName.Gaussian, dict(guess=True, prefix="p_")),
               (ModelName.Linear, dict(slope=0.0, intercept=10.0, prefix="bg_")),
           ],
           "pol": [(ModelName.Cosine, dict(guess=True))],
       },
       averages=["temp"],
   )

   results = GroupFit(scans, params).fit()  # scans: Sequence[Scan], fit in this order

   for result in results:
       print(
           result.uuid.value,
           result.fits["detector"]["p_"]["center"],  # peak center
           result.fits["pol"].frequency,             # single-component fits expose parameters directly
           result.averages["temp"].mean,
           result.averages["temp"].std,
       )


Parameters
==========

``GroupFitParams`` says what to do with each scan's columns:

``x_col``
   The independent variable every fitted column is fit against. ``None`` uses each scan's
   ``tavimeta.default_axis[0]``.

``fits``
   ``{y column: model_dict}``. A peak and a sinusoid differ only in their model
   (``ModelName.Gaussian`` vs ``ModelName.Cosine``), so both go here. Each point is weighted by its
   counting error, ``sqrt(|y|)``.

``averages``
   Columns reduced to their per-scan mean and standard deviation (the spread during the scan)
   instead of being fit, e.g. sample temperature.

``chain`` (default ``True``)
   Start each scan after the first from the previous scan's fitted values instead of the group
   model's starting values or guess. See `Chained fits`_.

``overrides``
   Per-scan changes to the fit inputs. See `Per-scan overrides`_.

Column names are looked up in ``scan.data``; a digit-leading name such as ``2theta`` also finds a
column stored as ``_2theta``. A missing column raises ``KeyError`` listing the scan's columns.


Results
=======

``fit()`` returns one ``GroupFitResult`` per scan, in scan order, and keeps the list on
``GroupFit.results``. Each holds:

- ``uuid`` — the source scan's uuid.
- ``fits`` — ``{y column: FitResult}``; read components by prefix, e.g. ``fits["detector"]["p_"]``.
- ``averages`` — ``{column: ColumnAverage(mean, std)}``.


How a scan's model is built
===========================

``GroupFit.model_for(uuid, y_col)`` returns the model a scan's column is fit with. Each component of
the group model is built in three layers, later layers winning:

.. mermaid::

   flowchart LR
       G["Group model<br/>fits[y_col]"] --> S["Chained seed<br/>previous scan's fitted values"]
       S --> O["Override<br/>this scan's patch"]
       O --> F["Fit.fit"]

Components are matched across layers by their ``prefix`` (``""`` when omitted), the same key
``Fit.fit`` uses to namespace parameters and ``FitResult`` uses to return them.


Chained fits
------------

With ``chain=True`` each component of scan *n* starts from scan *n - 1*'s fitted values for the
same column and prefix. This tracks a peak that drifts or fades across a sweep, where every scan
started from the same guess would lose it.

- Only ``Fit.fit`` starting values are seeded: ``center``/``sigma``/``amplitude`` (peaks),
  ``slope``/``intercept`` (Linear), ``amplitude``/``frequency``/``phase``/``offset`` (Cosine).
  Derived values (``fwhm``, ``height``) and Voigt's ``gamma`` are not.
- A seeded component has ``guess`` dropped, since ``Fit.fit`` ignores explicit starting values on a
  guessed component.
- If the previous fit lacks that component, or any seed value is non-finite, the component keeps
  its group starting values.
- A seed negligible next to its parameter's scale in the scan's data (under 1e-6 of it, e.g. a
  center of 2e-13 on an axis spanning 5) is seeded as exactly 0. lmfit's finite-difference step is
  proportional to the starting value, so a near-zero seed would leave that parameter frozen.
- ``set`` options (fixed, bounded or tied parameters) still apply on top of the seed.
- The first scan always uses the group model.


Per-scan overrides
------------------

A scan that needs hand-tuning (a weak peak whose width should be fixed, a peak that needs bounds)
gets a **patch** for one component, keyed by scan uuid, y column and component prefix. Patches can be
given up front:

.. code-block:: python

   params = GroupFitParams(
       fits={"detector": [...]},
       overrides={
           UUID(value="s-30"): {
               "detector": {
                   "p_": {"set": {"sigma": {"value": 0.2, "vary": False}, "center": {"min": 0.4, "max": 0.8}}},
               },
           },
       },
   )

or added later, then applied by fitting again:

.. code-block:: python

   group = GroupFit(scans, params)
   group.fit()
   group.override(UUID(value="s-30"), "detector", "p_", {"set": {"sigma": {"value": 0.2, "vary": False}}})
   group.fit()

Merge rules:

- ``set`` merges **per parameter**: patching ``{"sigma": {"vary": False}}`` onto a group
  ``{"sigma": {"min": 0}}`` gives ``{"sigma": {"min": 0, "vary": False}}``.
- Every other key (``guess``, ``center``, ``slope``, ...) replaces the value from the group model or
  the chained seed.
- Successive ``override()`` calls on the same component accumulate by the same rules.
- On a guessed component, change a starting value with ``set={"center": {"value": ...}}``; a bare
  ``center`` would be ignored by ``Fit.fit``.

Up-front overrides go through ``override()`` when the ``GroupFit`` is built, so a scan uuid, column
or prefix the group does not have raises ``KeyError`` immediately. ``GroupFit`` works on a private
copy of the ``GroupFitParams`` it was given; ``override()`` never changes the caller's object.

Because a refit of one scan does not update the scans chained after it, there is no single-scan
refit: call ``fit()`` again after changing overrides.


GroupFitModel
=============

``tavi.backend.model.group_fit_model.GroupFitModel`` runs a group fit from the Group Fit window's
own state: each scan's fitting-panel ``FitSpec`` plus one shared y-vs-x column pair, packed in a
``GroupFitRequest``. It parses each spec exactly as a single fit does (range, blank fields guessed
from the data, FWHM to sigma, fixed and bounded parameters) by reusing ``FitModel``'s translation,
then hands the group to ``GroupFit``:

- The first scan's spec is the group model.
- With ``chain=True``, a later scan gets an override only for the parameters the user changed
  relative to the first scan, so the chained seed supplies the rest. A changed parameter whose
  value was left blank overrides only its fixed/bounds options.
- With ``chain=False``, every later scan overrides every component, so it is fit from its own spec.
- Each scan's x and y are trimmed to that scan's own fitting range before ``GroupFit`` sees them.
  A blank bound leaves that side open, so a blank range is the whole scan. ``averages`` still
  cover the whole scan.
- A scan that can't be fit (not loaded, bad range or column, or a background and peaks unlike the
  first scan's) is reported and left out, and listed without a result. Every problem in one run
  is gathered into a single error report.

The result is published as one ``SyncGroupFitEvent`` with a ``GroupFitOutcome`` per scan: the
spec, the ``FitResultSummary``, the curve, the fitted points and the column averages. This is
display only: the fit itself is not saved into the project.

``summarize_group_fit`` saves one fitted parameter against a column's per-scan average. It takes
a ``SummarizeGroupFitRequest``: a key of ``FitResultSummary.parameters()`` (e.g.
``"peak1 center"``), a column name, and each scan's ``FitResultSummary`` from the last group fit.
For every scan it pairs the parameter and its 1-sigma with the column's mean and standard
deviation over the whole scan, and publishes one ``SaveFitSummaryEvent`` carrying a
``FitSummary``. ``TaviProjectModel`` stores it in ``TaviData.fits`` and lists it under Fits. A
scan that is no longer loaded, lacks the column, or lacks the parameter is reported and left out.

``GroupFitPresenter`` connects it to the Group Fit window; see :doc:`../design/frontend/group_fit_flow`.


API Reference
=============

.. autoclass:: tavi.library.fit.group_fit.GroupFitParams
   :members:

.. autoclass:: tavi.library.fit.group_fit.GroupFit
   :members:

.. autoclass:: tavi.library.fit.group_fit.GroupFitResult
   :members:

.. autoclass:: tavi.library.fit.group_fit.ColumnAverage
   :members:
