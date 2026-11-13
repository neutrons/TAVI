Histogram
=========

.. class:: Histogram(tavi_data: TaviData, uuids: Sequence[UUID], columns: Sequence[str], bins: Any = None, range: Optional[Sequence[tuple[float, float]]] = None, multiplication: float = 1.0, tol: Sequence[float] = (0.01, 0.01))

   Append several raw scans, then bin them on ``columns = [x, y, z, norm]`` into
   a regular 2D intensity map. Lives in
   ``tavi.library.data.algorithms.histogram``; the result is a
   :class:`HistogramData`, defined in ``tavi.library.data.scan``.

   It is the 2D counterpart of ``EqualStepRebinOps`` (see :doc:`process_data`):
   the same append-then-bin approach and the same per-bin arithmetic, with the
   binning itself done by :func:`numpy.histogramdd`.

   :param tavi_data: The pool the origin uuids are looked up in.
   :type tavi_data: TaviData
   :param uuids: Origin scan uuids, at least one, each of which must be a ``RawScan``.
   :type uuids: Sequence[UUID]
   :param columns: Exactly four distinct column names, in the order
      ``[x, y, z, norm]`` — the two axes, the raw counts, and the monitor or
      time they were counted against.
   :type columns: Sequence[str]
   :param bins: Passed to ``numpy.histogramdd``: an int, ``[nx, ny]``, or
      ``[x_edges, y_edges]`` (forms may be mixed per axis). ``None`` reads a bin
      count per axis off the data — see `Choosing the bins`_.
   :type bins: Any
   :param range: ``[(xmin, xmax), (ymin, ymax)]``, passed to
      ``numpy.histogramdd``; points outside it are dropped. ``None`` reads the
      extent off the data.
   :type range: Optional[Sequence[tuple[float, float]]]
   :param multiplication: Scale applied to intensity and error, so that e.g.
      ``norm="monitor", multiplication=10000`` reads as counts per 10000 monitor
      counts.
   :type multiplication: float
   :param tol: ``[x_tol, y_tol]``, the tolerance within which values along each
      axis count as the same position when ``bins=None``. Ignored once ``bins``
      is given.
   :type tol: Sequence[float]

Overview
--------

.. code-block:: python

    from tavi.library.data.algorithms.histogram import Histogram

    hist = Histogram(tavi_data, [scan_175.uuid, scan_176.uuid], ["l", "e", "detector", "mcu"]).exec()

    hist.tavimeta.default_axis      # ("l", "e")
    hist.x_centers                  # [1.6499, 1.8001]  - one bin per constant-Q scan
    len(hist.y_centers)             # 21                - one bin per energy step
    np.shape(hist.intensity)        # (2, 21)

``exec()`` runs ``AppendOp`` over the origins, takes the four columns as float
arrays, computes each point's error as ``sqrt(z)``, and bins them. The result is
a fresh scan in its own right: it keeps the append's new uuid, its
``prov.contributing_scans`` and its ``+``-joined friendly name (suffixed
``_hist_norm<column>``, plus ``x<multiplication>`` when that is not 1). Like
``CombineOps``, ``exec()`` then **registers** it in ``TaviData.processed_scans``
through ``register()`` — kept separate so a map rebuilt elsewhere can be stored
without re-binning — and returns it as well.

.. mermaid::

    flowchart LR
        A["RawScans<br/>(uuids)"] --> B["AppendOp<br/>[x, y, z, norm]"]
        B --> C["x, y, z, norm arrays<br/>err = sqrt(z)<br/>drop non-finite rows"]
        C --> D{"bins?"}
        D -- None --> E["count positions<br/>per axis within tol"]
        E --> F
        D -- "int / [nx, ny]" --> F["pad range half a bin<br/>past extreme points"]
        D -- "edges or range given" --> G
        F --> G["np.histogramdd x4<br/>weights: z, err², norm, none"]
        G --> H["HistogramData"]

Per-bin arithmetic
------------------

Each bin combines its points the way one longer measurement would, exactly as
the 1D rebins do. For the points :math:`i` falling in a bin:

.. math::

   I = m \, \frac{\sum_i z_i}{\sum_i n_i}, \qquad
   \sigma = m \, \frac{\sqrt{\sum_i z_i}}{\sum_i n_i}

where :math:`n` is the norm column and :math:`m` is ``multiplication``. The
error is propagated as :math:`\sqrt{\sum \sigma_i^2}` with
:math:`\sigma_i = \sqrt{z_i}`, the Poisson error of raw counts, which reduces to
the expression above.

All four sums come from ``numpy.histogramdd`` on the same ``(N, 2)`` sample of
``(x, y)`` points, one call per weight:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - ``weights``
     - Per-bin result
   * - ``z``
     - :math:`\sum z` — summed counts
   * - ``err**2``
     - :math:`\sum \sigma_i^2`
   * - ``norm``
     - :math:`\sum n` — stored in ``tavimeta.normalization``
   * - none
     - number of points — stored as ``counts``

The edges returned by the first call are passed as ``bins`` to the other three,
so all four land on an identical grid whatever form ``bins`` and ``range`` were
given in.

Choosing the bins
-----------------

The map is always an equal-step grid laid by ``numpy.histogramdd``; ``bins`` and
``range`` only decide how it is laid.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Given
     - Grid
   * - ``bins=None`` (default)
     - Along each axis, values within ``tol`` of the first value of their group
       count as one position (the ``ToleranceRebinOps`` rule), and the axis gets
       as many bins as it has positions. Then as for a count below.
   * - a count, no ``range``
     - The extent is widened half a bin past the extreme points, so the
       outermost bins are **centred** on them.
   * - a count and a ``range``
     - ``numpy.histogramdd`` as is: ``range`` gives the outer edges.
   * - edge arrays
     - Used as given; ``range`` is ignored for that axis.

Points outside the grid are dropped. A point exactly on the last edge falls in
the last bin, as in ``numpy.histogramdd`` and the 1D equal-step rebin.

For two constant-Q energy scans at ``l = 1.65`` and ``1.8``:

.. code-block:: python

    Histogram(..., bins=[2, 21])                          # x centres at the extreme readbacks, 1.6499 / 1.8001
    Histogram(..., bins=[np.array([1.575, 1.725, 1.875]), 21])   # x centres exactly 1.65 / 1.8

HistogramData
-------------

.. class:: HistogramData

   A ``Scan`` (``BaseModel, Scan``, every field writable, like
   ``ProcessedScan``) carrying a 2D intensity map. Every grid is shaped
   ``(n_x, n_y)``, ``numpy.histogramdd``'s own order.

   .. list-table::
      :header-rows: 1
      :widths: 30 25 45

      * - Field
        - Type
        - Holds
      * - ``uuid``, ``data``, ``metadata``, ``prov``
        - as ``Scan``
        - From the append: fresh uuid, the appended points, origins in provenance.
      * - ``tavimeta.default_axis``
        - ``(str, str)``
        - The x and y column names.
      * - ``tavimeta.normalization``
        - ``(str, np.ndarray)``
        - The norm column and each bin's summed norm, an ``(n_x, n_y)`` array.
      * - ``intensity_name``
        - ``str``
        - The z column the intensity was binned from.
      * - ``x_edges``, ``y_edges``
        - ``np.ndarray``
        - ``n + 1`` bin edges per axis.
      * - ``x_centers``, ``y_centers``
        - ``np.ndarray``
        - Bin midpoints, ``(edges[:-1] + edges[1:]) / 2``.
      * - ``intensity``, ``error``
        - ``np.ndarray``
        - Normalized intensity and its error; ``NaN`` where no point fell.
      * - ``counts``
        - ``np.ndarray`` (int)
        - Number of points in each bin.

   ``TaviMetadata.normalization`` is typed
   ``Optional[Tuple[str, float] | Tuple[str, np.ndarray]]`` to admit the grid
   (``TaviMetadata`` allows arbitrary types for it); 1D scans keep storing
   ``(column, factor)``.

Plotting
~~~~~~~~

The grids are numpy arrays; transpose for ``pcolormesh``, which wants rows
along y:

.. code-block:: python

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(hist.x_edges, hist.y_edges, hist.intensity.T, shading="flat")
    fig.colorbar(mesh, ax=ax, label=f"{hist.intensity_name} / {hist.tavimeta.normalization[0]}")
    ax.set_xlabel(hist.tavimeta.default_axis[0])
    ax.set_ylabel(hist.tavimeta.default_axis[1])

``NaN`` bins are left blank. A single cut is a row of the grid, e.g.
``ax.errorbar(hist.y_centers, hist.intensity[0], yerr=hist.error[0])``.

Validation
----------

``validate()`` (called by ``exec()``) raises:

- ``ValueError`` if ``columns`` is not exactly four distinct names; if
  ``multiplication`` is not positive and finite; if ``tol`` is not a pair of
  non-negative finite numbers; if no uuids were given; or if an origin is not a
  ``RawScan``.
- ``KeyError`` if an origin does not carry one of the columns (from ``AppendOp``).

``exec()`` additionally raises ``ValueError`` if no row has finite values in all
four columns, or if ``bins=None`` and ``range`` holds none of the points. Rows
with a non-finite value anywhere are dropped rather than binned.

Key Design Decisions
--------------------

Weights carry the physics, ``sample`` only positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``numpy.histogramdd`` bins along every column of ``sample``, so only ``x`` and
``y`` go there; adding ``z`` would make a 3D histogram of point counts. Counts,
squared errors and norm are summed into the bins through ``weights``, and the
intensity is formed from those sums afterwards. Summing counts and norm
separately — rather than averaging per-point ratios — is what makes a bin
behave like one longer measurement, matching the 1D rebins.

Raw scans only, error from counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The error is not a column the caller passes; it is ``sqrt(z)`` for every point.
That is the right error only for raw counts, so — like ``EqualStepRebinOps`` and
``ToleranceRebinOps`` — every origin must be a ``RawScan``.

Empty bins are ``NaN``
~~~~~~~~~~~~~~~~~~~~~~

A bin no point fell in reports ``NaN`` intensity and error, not 0. Zero is a
measurement; an empty bin is an absence. matplotlib leaves ``NaN`` cells blank.

Outer bins are centred on the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Left to itself, ``numpy.histogramdd`` puts the outer *edges* at the extreme
data values. Scans at ``l = 1.65`` and ``1.8`` binned into two would then sit on
the grid's edges, with centres at ``1.6875`` and ``1.7625``. For a count with no
``range``, the extent is widened half a bin either side so the extreme points
are bin centres — the same fallback ``EqualStepRebinOps`` uses when ``start`` or
``finish`` is left out. Explicit edges or a ``range`` are honoured as given.

``bins=None`` counts positions, it does not invent a grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Motor readback jitters (``1.6499``, ``1.65``, ``1.6501``), so matching values
exactly would split one setpoint into several narrow bins. ``bins=None``
therefore groups values within ``tol`` and only *counts* the groups, handing
that count to the ordinary equal-step path. The map stays a regular grid that
``pcolormesh`` can draw; the cost is that unevenly spaced positions do not each
land on a bin centre — pass edges for that.

Grids are numpy arrays
~~~~~~~~~~~~~~~~~~~~~~

``HistogramData`` keeps every grid — and the summed-norm grid in
``tavimeta.normalization`` — as the numpy array ``numpy.histogramdd`` produced,
rather than converting to nested lists as ``ScanData`` stores its columns. Each
grid's ``(n_x, n_y)`` shape travels with it and plotting needs no conversion.
``HistogramData`` and ``TaviMetadata`` allow arbitrary types to admit them.

The scan still deep-copies, which every ``EventBroker`` publish does, but it no
longer serializes to JSON: once a map is stored in a project, saving it will
need the arrays converted on the way out and back on load.

Stored beside the 1D combinations, run outside ``CombineOps``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Histogram`` is not a ``ProcessOps`` subclass, whose ``exec`` contract returns
a ``ProcessedScan``, and ``CombineOps`` has no 2D mode. The Combine panel's
Histogram tab reaches it through ``CombineModel`` instead, which turns each
axis's Start/Stop/Step into explicit bin edges (a blank Step gives one bin per
measured position, values within ``POSITION_TOL`` counting as one, with edges
midway between neighbours) and runs it; ``Histogram.exec()`` stores the result
in ``TaviData.processed_scans`` — typed ``ProcessedScan | HistogramData`` — so
a map is listed, renamed and focused like any other combination. Clicking it
publishes ``FocusHistogramEvent`` and the plotter draws it with
``pcolormesh`` (see :doc:`frontend/visualization_flow`).
