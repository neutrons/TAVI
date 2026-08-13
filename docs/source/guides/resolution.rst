.. _resolution:

++++++++++
Resolution
++++++++++

The resolution stack computes the 4D (:math:`Q_x, Q_y, Q_z, E`) resolution
matrix of a triple-axis spectrometer at a given :math:`(hkl, E_i, E_f)`, and
projects it into the frame the user wants to look at.

Three classes are involved:

* ``Resolution`` -- manager. Selects the model, assembles the instrument/sample
  geometry, and returns the resolution matrix together with the normalisation
  factor ``r0``.
* ``CooperNathans`` -- the model itself. Produces the matrix in the *local*
  :math:`Q` frame.
* ``ResolutionEllipsoid`` -- projects the local matrix into an ``hkle`` frame and
  extracts coherent / incoherent FWHM along a chosen axis.

The high-level entry points are :ref:`TAS.resolution_bar and TAS.browse
<triple_axis>`, which drive all three.


Overview
========

Source modules:

.. code-block:: text

   tavi.library.resolution.resolution
   tavi.library.resolution.ellipsoid
   tavi.library.resolution.cooper_nathans

Available models:

.. autoclass:: tavi.library.resolution.resolution.ResolutionModel
   :members:
   :undoc-members:

Manager:

.. autoclass:: tavi.library.resolution.resolution.Resolution
   :members:
   :undoc-members:

Ellipsoid and projections:

.. autoclass:: tavi.library.resolution.ellipsoid.ResolutionEllipsoid
   :members:
   :undoc-members:


.. _resolution_frame:

Resolution Frames
=================

The frame is selected with the ``resolution_frame`` argument, accepted by
``Resolution``, ``ResolutionEllipsoid``, ``TAS.resolution_bar`` and
``TAS.browse``. Two kinds of value are supported.

Local :math:`Q` frame
---------------------

.. code-block:: python

   resolution_frame = "local"

The matrix is returned exactly as the model produces it, in the frame attached
to the scattering triangle. No projection and no rotation matrix are needed, so
this is the cheapest option and the default for ``TAS.browse`` and
``TAS.resolution_bar``.

Axis order:

=====  ==================================================================
Index  Meaning
=====  ==================================================================
0      :math:`Q` parallel (along :math:`\vec Q`), in :math:`\AA^{-1}`
1      :math:`Q` perpendicular, in plane, in :math:`\AA^{-1}`
2      :math:`Q` vertical, in :math:`\AA^{-1}`
3      Energy transfer, in meV
=====  ==================================================================

``resolution_frame=None`` is treated the same as ``"local"``.

``hkle`` frames
---------------

.. code-block:: python

   resolution_frame = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e")   # plain HKL
   resolution_frame = ((1, 1, 0), (1, -1, 0), (0, 0, 1), "e")  # rotated frame
   resolution_frame = ((1, 0, 0), (0, 0, 1), (0, 1, 0), "e")   # reordered axes

A 4-tuple of three reciprocal-lattice vectors :math:`(w_1, w_2, w_3)` plus the
literal ``"e"`` for the energy axis. The local matrix is transformed with

.. math::

   M_{hkle} = C^T M_{local} C, \qquad
   C_{3\times3} = 2\pi \, R_{lab \leftarrow local} \, R \, UB \, W

where :math:`R` is the sample rotation matrix, :math:`UB` comes from the
sample's oriented lattice, and :math:`W = [w_1\ w_2\ w_3]` is the column matrix
of the projection vectors. The plain HKL frame is the special case
:math:`W = I` and is short-circuited in ``project_to_frame``.

Axis order is the order of the tuple: index 0, 1, 2 are :math:`w_1, w_2, w_3`
(in r.l.u.) and index 3 is energy in meV. Reordering the vectors reorders the
axes -- ``((1, 0, 0), (0, 0, 1), (0, 1, 0), "e")`` puts :math:`(0, 0, 1)` on
index 1 and :math:`(0, 1, 0)` on index 2.

The three vectors must be non-coplanar; ``get_axes_angles`` raises
``ValueError`` otherwise. The angles between the basis vectors are computed from
the reciprocal lattice, so a non-orthogonal frame is drawn on a skewed grid by
``PlotResolution``.


.. _projection_axis:

Projection Axis
===============

``projection_axis`` selects which axis of the (possibly projected) matrix the
FWHM is taken along. It indexes the frame chosen above, so its meaning depends
on ``resolution_frame``:

* ``resolution_frame="local"``, ``projection_axis=0`` -- FWHM along
  :math:`\vec Q`; ``1`` -- in-plane perpendicular; ``3`` -- energy.
* ``resolution_frame=((1, 0, 0), (0, 0, 1), (0, 1, 0), "e")``,
  ``projection_axis=3`` -- energy FWHM in that frame; ``0`` -- FWHM along
  :math:`(1, 0, 0)` in r.l.u.

Two widths are available:

``coh_fwhm(projection_axis)``
   Coherent FWHM -- a *cut*. Read straight off the diagonal element,
   :math:`\mathrm{FWHM} = 2\sqrt{2\ln 2} / \sqrt{M_{ii}}`. This is what the
   resolution bar draws.

``incoh_fwhm(projection_axis)``
   Incoherent FWHM -- an *integral*. The other three axes are removed with
   ``quadric_proj`` before the width is taken, giving the broader,
   integrated width. ``TAS.resolution_bar`` computes it for ``hkle`` frames and
   currently prints it; it is not plotted yet.

For 2D slices use ``Resolution.get_ellipse(res_mat, ellipse_axes=(0, 1))``,
which keeps two axes and returns the 2x2 matrix plus the angle between the two
kept basis vectors. ``Resolution.get_projection(res_mat, projection_axes)``
generalises this to any subset of axes -- pass ``projection_axes`` as a tuple of
the indices to keep. In both, ``PROJECTION=False`` (the default) takes a cut by
deleting rows/columns, while ``PROJECTION=True`` integrates via
``quadric_proj``.


Return Value
============

``Resolution.get_resolution(hkl, ei, ef, rot_mat=None)`` always returns a
2-tuple:

.. code-block:: python

   res_mat, r0 = resolution.get_resolution(hkl=(0, 0, 2), ei=14.5, ef=14.5)

``res_mat`` is the 4x4 matrix in the requested frame and ``r0`` the
normalisation factor, kept alongside the matrix so intensities can be
normalised downstream.


.. _resolution_rotation_matrix:

Rotation Matrix
===============

Projecting into an ``hkle`` frame needs the sample rotation matrix :math:`R`.
There are two ways to get it, and which one is used depends on the caller:

``Goniometer.r_mat(angles)`` (explicit)
   Built from the measured motor angles. ``TAS.resolution_bar`` uses this for
   ``hkle`` frames, passing ``rot_mat`` into ``get_resolution``. It requires

   * the goniometer ``type`` in the instrument JSON to match one of the
     implemented cases exactly -- ``"Y,-Z,X"`` (Huber table, uses ``omega``,
     ``sgl``, ``sgu``) or ``"Y,Z,Y,bisect"`` (four-circle, uses ``omega``,
     ``chi``, ``phi``). The strings carry **no spaces**. ``hb1.json``,
     ``hb1a_tas.json``, ``hb3.json`` and ``cg4c.json`` declare ``"Y,-Z,X"``;
     ``hb1a_4c.json`` declares ``"Y,Z,Y,bisect"``. Note that
     ``Goniometer.__init__`` defaults to ``"Y, -Z, X"`` *with* spaces, which
     matches neither -- a hand-built goniometer must pass ``type`` explicitly.
   * the ``MotorAngles`` of the data point to carry the keys the goniometer
     expects. ``ORNLSpiceLoader.get_data_point_closest_to_center`` populates
     ``two_theta``, ``omega``, ``sgl`` and ``sgu`` for SPICE scans that report
     ``s1``/``s2``.

   Anything else hits ``case _`` in ``Goniometer.r_mat`` and raises
   ``ValueError("Unknown mode.")``, which ``TAS.resolution_bar`` re-raises as
   ``ValueError("Rotation matrix can not be calculated")``.

``OrientedLattice.rot_matrix_with_minimal_tilt`` (implicit)
   Used when ``rot_mat=None`` is passed to ``get_resolution`` for an ``hkle``
   frame. Derives the orientation from the UB matrix with the smallest tilt that
   reaches ``hkl``. Needs ``UB``, ``plane_normal`` and ``in_plane_ref`` to be set
   on ``sample.ol``; a ``ValueError`` pointing at those fields is raised
   otherwise.

The local frame needs neither, so ``rot_mat`` is ignored when
``resolution_frame="local"``.


Minimal Example
===============

Resolution at a single point, in a rotated ``hkle`` frame:

.. code-block:: python

   from tavi.library.Instrument.instrument import Instrument, InstrumentCatalog
   from tavi.library.experiment.experiment import Experiment
   from tavi.library.experiment.enum import FixedEnergyMode
   from tavi.library.resolution.ellipsoid import ResolutionEllipsoid
   from tavi.library.resolution.resolution import Resolution, ResolutionModel

   experiment = Experiment(mode=FixedEnergyMode.FIX_Ef, fixed_energy=14.503)
   experiment.load_folder("/path/to/exp1091/Datafiles")
   sample = experiment.create_sample("/path/to/UBConf/UBmatrix.dat")
   instrument = Instrument(instrument_catalog=InstrumentCatalog.HB1A4C)

   frame = ((1, 0, 0), (0, 0, 1), (0, 1, 0), "e")
   resolution = Resolution(
       model=ResolutionModel.CooperNathans,
       instrument=instrument,
       sample=sample,
       experiment=experiment,
       resolution_frame=frame,
   )

   res_mat, r0 = resolution.get_resolution(hkl=(0, 0, 2), ei=14.503, ef=14.503)

   # Energy FWHM (index 3) in this frame.
   ellipsoid = ResolutionEllipsoid(res_mat, resolution_frame=frame)
   print(ellipsoid.coh_fwhm(projection_axis=3), ellipsoid.incoh_fwhm(projection_axis=3))

   # 2D slice in the (w1, energy) plane, plus the angle between those axes.
   res_2d, axes_angle = resolution.get_ellipse(res_mat, ellipse_axes=(0, 3))
