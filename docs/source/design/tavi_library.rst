Triple-Axis Library
===================

The ``tavi.library`` triple-axis stack provides the in-memory representation of a
triple-axis spectrometer (TAS) and the resolution calculations performed on top
of it. This page documents the modules added in the Cooper-Nathans resolution
branch.

Core Concepts
-------------

The stack is organised around five classes:

- **Instrument** (``Instrument``)
  Container holding the spectrometer components (monochromator, analyzer,
  collimators, goniometer).

- **Components** (``Crystal``, ``Collimators``, ``Goniometer``)
  Reusable hardware abstractions used by ``Instrument`` and consumed by the
  resolution calculation.

- **Experiment** (``Experiment``)
  Wraps a ``Scan``, tracks the fixed-energy mode (``FixedEnergyMode``), and
  exposes the geometric quantities derived from the data (two-theta, psi,
  peak indices).

- **Resolution** (``Resolution``, ``CooperNathans``, ``ResolutionEllipsoid``)
  Computes the 4D resolution matrix and ``r0`` normalisation factor and
  projects them into user-selected frames.

- **TAS** (``TAS``)
  Top-level façade tying ``Instrument``, ``Sample``, and ``Experiment``
  together for downstream calculations.

Architecture Overview
---------------------

Components by responsibility:

- ``tavi.library.Instrument.instrument.Instrument`` -- aggregate spectrometer
- ``tavi.library.component.crystal.Crystal`` -- monochromator/analyzer crystals
- ``tavi.library.component.collimators.Collimators`` -- soller-slit divergences
- ``tavi.library.component.goniometer.Goniometer`` -- sample stage
- ``tavi.library.experiment.experiment.Experiment`` -- scan-derived geometry
- ``tavi.library.resolution.resolution.Resolution`` -- resolution manager
- ``tavi.library.resolution.cooper_nathans.CooperNathans`` -- Cooper-Nathans model
- ``tavi.library.resolution.ellipsoid.ResolutionEllipsoid`` -- 4D ellipsoid + projection
- ``tavi.library.tas.triple_axis.TAS`` -- top-level entry point

Instrument
----------

``Instrument`` aggregates the optical and mechanical components of the
spectrometer.

Constructor arguments:

- ``monochromator`` (``Crystal``) -- pre-sample crystal
- ``analyzer`` (``Crystal``) -- post-sample crystal
- ``collimators`` (``Collimators``) -- four sets of horizontal and vertical divergences
- ``goniometer`` (``Goniometer``) -- sample stage with scattering sense

Example:

.. code-block:: python

    from tavi.library.Instrument.instrument import Instrument
    from tavi.library.component.crystal import Crystal
    from tavi.library.component.collimators import Collimators
    from tavi.library.component.goniometer import Goniometer

    instrument = Instrument(
        monochromator=Crystal(crystal="PG002", sense="+"),
        analyzer=Crystal(crystal="PG002", sense="-"),
        collimators=Collimators(pre_mono_h=40, pre_sample_h=100),
        goniometer=Goniometer(s2_sense="+"),
    )

Instrument parameter files (``hb1.json``, ``hb1a.json``, ``hb3.json``,
``cg4c.json``) live in ``tavi/library/Instrument/instrument_params/`` and
provide default configurations for the HFIR instruments.

Components
----------

Crystal
~~~~~~~

``Crystal`` represents a monochromator or analyzer crystal. It carries a
mosaic (horizontal/vertical), a crystal type, and a scattering sense.

Constructor arguments:

- ``mosaic_h``, ``mosaic_v`` (``float``) -- mosaic FWHM in arc-minutes
- ``crystal`` (``str``) -- e.g. ``"PG002"``, ``"Cu220"``, ``"Heusler"``
- ``sense`` (``"+"`` or ``"-"``) -- scattering sense

The supported crystals and their d-spacings (in angstrom) are listed in
``crystal_d`` at the top of ``crystal.py`` (``PG002``, ``PG004``, ``Cu111``,
``Cu220``, ``Ge111``, ``Ge220``, ``Ge311``, ``Ge331``, ``Be002``, ``Be110``,
``Heusler``). If the crystal is not listed here, a custom string can be set
for crystal and specific d-spacing can be manually put in.

Collimators
~~~~~~~~~~~

``Collimators`` holds the horizontal and vertical divergence (in arc-minutes,
stored internally and exposed as radians) for the four soller slits:

- ``pre_mono_h`` / ``pre_mono_v``
- ``pre_sample_h`` / ``pre_sample_v``
- ``post_sample_h`` / ``post_sample_v``
- ``post_ana_h`` / ``post_ana_v``

Properties return values in radians (divided by 60 and converted via
``np.radians``) for direct use in matrix calculations.

Goniometer
~~~~~~~~~~

The goniometer carries the scattering sense (``+`` or ``-``) used to sign
``two_theta`` and ``psi`` in geometric calculations.

Experiment
----------

``Experiment`` wraps an optional ``Scan``, records which energy is held fixed
during the scan, and exposes the geometric quantities needed by the
resolution machinery.

Constructor arguments:

- ``scan`` (``Scan``) -- optional scan payload
- ``mode`` (``FixedEnergyMode``) -- ``FixedEnergyMode.FIX_Ei`` or
  ``FixedEnergyMode.FIX_Ef``; defaults to ``FixedEnergyMode.FIX_Ef``
- ``fixed_energy`` (``float``) -- value of the fixed energy in meV

``FixedEnergyMode`` is an ``Enum`` with exactly two members; only enum
members are accepted (raw strings are not).

Methods:

.. code-block:: python

    def get_two_theta(q_norm: float, ei: float, ef: float) -> float
        """Scattering angle from ki, kf, and |Q|."""

    def get_psi(q_norm: float, ei: float, ef: float) -> float
        """Angle between ki and Q (sign opposite of s2)."""

    def set_fixed_energy(e: float) -> None
        """Assign ``e`` to ``ei`` or ``ef`` based on ``mode``."""

    def get_ei_ef(e: float) -> tuple[float, float]
        """Return ``(ei, ef)`` given the complementary energy ``e``."""

Both angles are derived from a triangle solve using ``ki = SE2K(ei)`` and
``kf = SE2K(ef)``. The scattering triangle in k-space (``Q = ki - kf``) is:

.. image:: ../images/scattering_triagnle.png
   :alt: Scattering triangle in k-space
   :align: center

    2θ : scattering angle, between ki and kf (returned by get_two_theta)
    ψ  : angle between ki and Q             (returned by get_psi)

Sign conventions are applied by the caller (``Resolution``) based on
instrument/goniometer sense.

Example:

.. code-block:: python

    from tavi.library.experiment.experiment import Experiment, FixedEnergyMode

    experiment = Experiment(mode=FixedEnergyMode.FIX_Ef, fixed_energy=4.8)

Resolution
----------

``Resolution`` is the entry point for resolution calculations. It selects a
model (currently ``"Cooper-Nathans"``), gathers geometry from the experiment
and instrument, and optionally projects the resulting matrix into user
selected axes.

Constructor arguments:

- ``model`` -- ``"Cooper-Nathans"``
- ``instrument`` (``Instrument``)
- ``sample`` (``Sample``)
- ``experiment`` (``Experiment``)
- ``scan_idx``, ``pt_idx`` (``int``) -- scan/point indices
- ``axes`` (``tuple``) -- projection axes; last entry is ``"e"`` for energy

Key methods:

.. code-block:: python

    def get_resolution(
        hkl: tuple[float, float, float],
        ei: float,
        ef: float,
        rot_mat: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, float]
        """Resolution matrix + r0 at ``hkl``, optionally projected via ``rot_mat``."""

    def get_ellipse(
        res_mat: np.ndarray,
        ellipse_axes: tuple[int, int] = (0, 1),
        PROJECTION: bool = False,
    ) -> np.ndarray
        """Reduce a 4D resolution matrix to a 2D ellipse by slice or projection."""

    def r_matrix_with_minimal_tilt(
        hkl: tuple[float, float, float],
        ei: float,
        ef: float,
    ) -> np.ndarray
        """Rotation matrix bringing the scattering plane in with minimal tilt."""

When ``axes`` is ``None``, ``get_resolution`` returns the raw matrix in the
local Q frame. Otherwise it constructs a ``ResolutionEllipsoid`` and calls
``project_to_frame`` to produce the matrix in the requested axes; in that
mode ``rot_mat`` is required.

CooperNathans
~~~~~~~~~~~~~

``CooperNathans`` implements the Popovici 1975 formulation. It builds the
matrices ``G`` (collimator divergence), ``F`` (mosaic), ``C`` (crystal
geometry), and ``A``/``B`` (lab-frame projection) used to assemble the 4D
resolution matrix and the ``r0`` normalisation factor.

Class attributes:

- ``NUM_MONO``, ``NUM_ANA`` -- one each
- ``NUM_COLLS = 4`` -- four soller slits
- ``IDX_COLLS`` -- index map for horizontal/vertical entries in ``G``

The main entry point is:

.. code-block:: python

    def resolution_matrix(
        instrument: Instrument,
        sample: Sample,
        q_norm: float,
        ki: float,
        kf: float,
        psi: float,
        two_theta: float,
        theta_m: float,
        theta_a: float,
    ) -> tuple[np.ndarray, float]

ResolutionEllipsoid
~~~~~~~~~~~~~~~~~~~

``ResolutionEllipsoid`` wraps the 4D resolution matrix and the ``r0`` normalisation
along with the projection axes. ``project_to_frame(r_mat, psi, ub)`` rotates
the matrix into the requested ``hkle`` frame; when the requested axes differ
from the default ``((1,0,0),(0,1,0),(0,0,1),"e")``, the matrix is further
transformed by the user-supplied basis ``W`` and the rows/columns are
permuted to place the ``"e"`` axis at the requested position.

TAS
---

``TAS`` is the top-level façade combining ``Instrument``, ``Sample``, and
``Experiment``. It exposes ``get_resolution_at_hkle`` for resolution lookups
and the R-matrix helpers used by ``Resolution``.

Example:

.. code-block:: python

    from tavi.library.tas.triple_axis import TAS

    tas = TAS(instrument=instrument, sample=sample, experiment=experiment)
    res_mat, r0 = tas.get_resolution_at_hkle(
        res_model="Cooper-Nathans",
        hkle=(0, 0, 3, 0.0),
        frame=None,
    )

Usage Examples
--------------

Resolution at a Single Point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from tavi.library.resolution.resolution import Resolution

    res = Resolution(
        model="Cooper-Nathans",
        instrument=instrument,
        sample=sample,
        experiment=experiment,
        axes=None,
    )
    res_mat, r0 = res.get_resolution(hkl=(0, 0, 3), ei=4.8, ef=4.8)

Projecting onto a User Frame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    res = Resolution(
        model="Cooper-Nathans",
        instrument=instrument,
        sample=sample,
        experiment=experiment,
        axes=((1, 1, 0), (0, 0, 1), (1, -1, 0), "e"),
    )
    rot_mat = res.r_matrix_with_minimal_tilt(hkl=(0, 0, 3), ei=4.8, ef=4.8)
    res_proj = res.get_resolution(hkl=(0, 0, 3), ei=4.8, ef=4.8, rot_mat=r_mat)

Reducing to a 2D Ellipse
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    res_2d_slice = res.get_ellipse(res_mat, ellipse_axes=(0, 3), PROJECTION=False)
    res_2d_proj  = res.get_ellipse(res_mat, ellipse_axes=(0, 3), PROJECTION=True)

``PROJECTION=False`` slices the matrix (integrates the unselected axes out by
deletion), while ``PROJECTION=True`` projects via the quadric-projection
identity (``quadric_proj``).

Plot Ellipse
------------

``tavi.library.plot.plot_ellipse`` provides two public objects for rendering
2D resolution ellipses on a skewed (non-orthogonal) grid:

- **``grid_helper(angle, nbins)``** -- builds a ``GridHelperCurveLinear`` that
  skews the Matplotlib axes by ``angle`` degrees, so that non-orthogonal
  reciprocal-space axes are drawn correctly.
- **``Plot``** -- accumulates one or more ``EllipseEntry`` objects and renders
  them together on a single axes.

``EllipseEntry`` dataclass fields:

- ``patch`` (``Ellipse``) -- the Matplotlib patch to draw
- ``x_extent``, ``y_extent`` (``float``) -- half-widths of the bounding box
  in data coordinates, used to auto-scale the axes limits
- ``origin`` (``tuple[float, float]``) -- centre of the ellipse in data coordinates

``Plot`` API
~~~~~~~~~~~~

.. code-block:: python

    class Plot:
        def __init__(self, axes_angle: float) -> None: ...
        """
        axes_angle: angle (degrees) between the two plot axes.
        """

        @staticmethod
        def create_ellipse(
            mat: np.ndarray,
            origin: tuple[float, float] = (0.0, 0.0),
            **kwargs,
        ) -> EllipseEntry: ...
        """
        Compute the FWHM ellipse from a 2x2 resolution matrix.
        Eigen-decomposition gives the semi-axes lengths and orientation angle.
        sig2fwhm is applied so dimensions are full-width values.
        """

        def add_ellipse(
            self,
            mat: np.ndarray,
            origin: tuple[float, float] = (0.0, 0.0),
            **kwargs,
        ) -> EllipseEntry: ...
        """Convenience wrapper: create_ellipse + append to the draw queue."""

        def plot(
            self,
            ax: Optional[Axes] = None,
            pad: float = 1.1,
            show: bool = True,
        ) -> Axes: ...
        """
        Draw all queued ellipses.  Creates a new figure with a skewed grid
        if ax is None.  Applies the shear transform so that patches are
        rendered in the tilted coordinate system.  A legend is shown
        automatically when any patch carries a visible label.
        """

Example
~~~~~~~

.. code-block:: python

    from tavi.library.plot.plot_ellipse import Plot

    # res_2d is a 2×2 slice/projection from res.get_ellipse(...)
    p = Plot(axes_angle=60.0)          # 60° between the two axes
    p.add_ellipse(res_2d, label="(0,0,3)")
    p.add_ellipse(res_2d_proj, origin=(0.05, 0.0), linestyle="--", label="projected")
    ax = p.plot(show=True)

Multiple ellipses (e.g. along a scan) can be added with different ``origin``
values before a single call to ``plot``.

Design Characteristics
----------------------

- **Composition over inheritance** -- ``Instrument`` aggregates components.
- **Model selection at runtime** -- ``Resolution`` dispatches on ``model``;
  Cooper-Nathans is currently the only implementation.
- **Frame-agnostic core** -- resolution matrix is computed in local Q;
  projection is a separable step via ``ResolutionEllipsoid``.
- **Sign conventions centralised** -- ``Resolution`` applies the goniometer/
  monochromator/analyzer ``sense`` to angles before handing off to the model.

Future Considerations
---------------------

- Additional resolution models (Popovici, Eckold-Sobolev, Monte Carlo).
- Loading instrument defaults from the ``instrument_params/*.json`` files.
- Vectorised resolution evaluation across a scan.
- Caching of intermediate matrices when ``hkl`` varies but ei/ef are fixed.
