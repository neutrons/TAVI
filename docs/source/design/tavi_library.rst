Triple-Axis Library
===================

The ``tavi.library`` triple-axis stack provides the in-memory representation of a
triple-axis spectrometer (TAS) and the resolution calculations performed on top
of it. This page documents the modules added in the Cooper-Nathans resolution
branch.

Core Concepts
-------------

The stack is organised around four concerns:

- **Instrument** (``Instrument``)
  Container holding the spectrometer components (monochromator, analyzer,
  collimators, goniometer) and the fixed-energy mode.

- **Components** (``Crystal``, ``Collimators``, ``Goniometer``)
  Reusable hardware abstractions used by ``Instrument`` and consumed by the
  resolution calculation.

- **Experiment** (``Experiment``)
  Wraps a ``Scan`` and exposes the geometric quantities derived from the data
  (two-theta, psi, peak indices).

- **Resolution** (``Resolution``, ``CooperNathans``, ``ResoEllipsoid``)
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
- ``tavi.library.resolution.cooper_nathan.CooperNathans`` -- Cooper-Nathans model
- ``tavi.library.resolution.ellipsoid.ResoEllipsoid`` -- 4D ellipsoid + projection
- ``tavi.library.tas.triple_axis.TAS`` -- top-level entry point

Instrument
----------

``Instrument`` aggregates the optical and mechanical components of the
spectrometer and tracks which energy is held fixed during a scan.

Constructor arguments:

- ``monochromator`` (``Crystal``) -- pre-sample crystal
- ``analyzer`` (``Crystal``) -- post-sample crystal
- ``collimators`` (``Collimators``) -- four soller-slit divergences
- ``goniometer`` (``Goniometer``) -- sample stage with scattering sense
- ``mode`` (``"fix_ei"`` or ``"fix_ef"``) -- which energy is held fixed
- ``ei_or_ef`` (``float``) -- the fixed energy value, in meV

Methods:

.. code-block:: python

    def set_ei_or_ef(e: float) -> None
        """Assign the fixed energy on the side selected by ``mode``."""

    def get_ei_ef(e: float) -> tuple[float, float]
        """Resolve the complementary energy from energy transfer ``e``."""

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
        mode="fix_ef",
        ei_or_ef=4.8,
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
``Heusler``).

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

``Experiment`` wraps an optional ``Scan`` and exposes the geometric quantities
needed by the resolution machinery.

Methods:

.. code-block:: python

    def get_two_theta(q_norm: float, ei: float, ef: float) -> float
        """Scattering angle from ki, kf, and |Q|."""

    def get_psi(q_norm: float, ei: float, ef: float) -> float
        """Angle between ki and Q (sign opposite of s2)."""

Both angles are derived from a triangle solve using ``ki = SE2K(ei)`` and
``kf = SE2K(ef)``. Sign conventions are applied by the caller (``Resolution``)
based on instrument/goniometer sense.

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
        r_mat: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, float]
        """Resolution matrix + r0 at ``hkl``, optionally projected via ``r_mat``."""

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
local Q frame. Otherwise it constructs a ``ResoEllipsoid`` and calls
``project_to_frame`` to produce the matrix in the requested axes; in that
mode ``r_mat`` is required.

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

ResoEllipsoid
~~~~~~~~~~~~~

``ResoEllipsoid`` wraps the 4D resolution matrix and the ``r0`` normalisation
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
    r_mat = res.r_matrix_with_minimal_tilt(hkl=(0, 0, 3), ei=4.8, ef=4.8)
    res_proj = res.get_resolution(hkl=(0, 0, 3), ei=4.8, ef=4.8, r_mat=r_mat)

Reducing to a 2D Ellipse
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    res_2d_slice = res.get_ellipse(res_mat, ellipse_axes=(0, 3), PROJECTION=False)
    res_2d_proj  = res.get_ellipse(res_mat, ellipse_axes=(0, 3), PROJECTION=True)

``PROJECTION=False`` slices the matrix (integrates the unselected axes out by
deletion), while ``PROJECTION=True`` projects via the quadric-projection
identity (``quadric_proj``).

Design Characteristics
----------------------

- **Composition over inheritance** -- ``Instrument`` aggregates components.
- **Model selection at runtime** -- ``Resolution`` dispatches on ``model``;
  Cooper-Nathans is currently the only implementation.
- **Frame-agnostic core** -- resolution matrix is computed in local Q;
  projection is a separable step via ``ResoEllipsoid``.
- **Sign conventions centralised** -- ``Resolution`` applies the goniometer/
  monochromator/analyzer ``sense`` to angles before handing off to the model.

Future Considerations
---------------------

- Additional resolution models (Popovici, Eckold-Sobolev, Monte Carlo).
- Loading instrument defaults from the ``instrument_params/*.json`` files.
- Vectorised resolution evaluation across a scan.
- Caching of intermediate matrices when ``hkl`` varies but ei/ef are fixed.
