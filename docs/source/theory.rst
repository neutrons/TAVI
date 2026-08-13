.. _theory:

Theoretical considerations
##########################

TAVI visualizes and analyses data from **triple-axis spectrometers** (TAS). This
page collects the conventions and relations the library implements, so the code
in ``tavi.library`` can be read against them.

Kinematics
==========

A triple-axis instrument selects an incident wavevector :math:`\vec k_i` with the
monochromator and a scattered wavevector :math:`\vec k_f` with the analyzer. The
momentum and energy transferred to the sample are

.. math::

    \begin{align}
        \vec Q &= \vec k_i - \vec k_f\\
        \Delta E &= E_i - E_f = \frac{\hbar^2}{2m_n}\left(k_i^2 - k_f^2\right)
    \end{align}

The conversion between energy in meV and wavevector in :math:`\AA^{-1}` is
implemented by :func:`tavi.library.experiment.utilities.SE2K`.

Most TAS experiments are run with one of the two energies held fixed. That choice
is modelled by :class:`tavi.library.experiment.enum.FixedEnergyMode`
(``FIX_Ei`` / ``FIX_Ef``); ``Experiment.get_ei_ef`` maps a measured energy
transfer back onto the :math:`(E_i, E_f)` pair using the mode's fixed energy.

Scattering triangle
===================

.. image:: images/scattering_triagnle.png
   :alt: Scattering triangle in k-space
   :align: center

The three vectors :math:`\vec k_i`, :math:`\vec k_f` and :math:`\vec Q` close a
triangle, so the in-plane geometry follows from the law of cosines. Two angles
are needed downstream:

.. math::

    \begin{align}
        Q^2 &= k_i^2 + k_f^2 - 2 k_i k_f \cos 2\theta\\
        \cos 2\theta &= \frac{k_i^2 + k_f^2 - Q^2}{2 k_i k_f}\\
        \cos \psi &= \frac{k_i^2 + Q^2 - k_f^2}{2 k_i Q}
    \end{align}

* :math:`2\theta` is the **scattering angle** between :math:`\vec k_i` and
  :math:`\vec k_f`, returned by ``Experiment.get_two_theta``.
* :math:`\psi` is the angle between :math:`\vec k_i` and :math:`\vec Q`,
  returned by ``Experiment.get_psi``.

Both are computed as bare magnitudes from the triangle. Sign is applied by the
caller: ``Resolution`` multiplies :math:`2\theta` by the goniometer's scattering
sense and :math:`\psi` by its negative, so that :math:`\psi` always carries the
opposite sign of ``s2``. Scattering can occur on either side of the incident
beam, and the goniometer's ``sense`` (``"+"`` / ``"-"``) is what fixes which.

Reciprocal space and the UB matrix
==================================

Sample orientation uses the standard Busing-Levy decomposition
:math:`UB = U \cdot B`, implemented in
:class:`tavi.library.geometry.oriented_lattice.OrientedLattice`:

* :math:`B` maps Miller indices to an orthonormal reciprocal-space frame attached
  to the crystal, built from the lattice parameters
  :math:`(a, b, c, \alpha, \beta, \gamma)`.
* :math:`U` is the orientation matrix rotating that crystal frame into the
  laboratory frame.

The magnitude of the momentum transfer at a reflection is then

.. math::

    |Q| = 2\pi\left|B\begin{pmatrix} h\\k\\l \end{pmatrix}\right|

(``OrientedLattice.q_norm_from_hkl``). The Cartesian convention is Mantid's /
the International Tables': :math:`Z` along the incoming beam, :math:`Y` up,
:math:`X` in the horizontal plane, right-handed
(https://docs.mantidproject.org/nightly/concepts/Lattice.html). SPICE uses a
different frame; :func:`tavi.library.experiment.utilities.spice_to_mantid` and
``mantid_to_spice`` convert between the two.

:math:`U` is determined from indexed peaks by
:class:`tavi.library.ubalgorithm.ub.UBAlgorithm`, which supports one peak plus a
scattering plane, two peaks, three peaks, or a least-squares fit over many peaks.
See :ref:`ub_algorithm`.

Instrument resolution
=====================

A triple-axis measurement does not sample a single point of
:math:`(\vec Q, E)` but a 4D ellipsoidal volume set by the collimations, the
monochromator/analyzer mosaics, and the scattering geometry. TAVI computes that
volume with the **Cooper-Nathans** formalism
(:class:`tavi.library.resolution.cooper_nathans.CooperNathans`), which returns a
:math:`4\times4` resolution matrix :math:`M` together with a normalisation
factor :math:`r_0`.

The matrix is produced in the *local* :math:`Q` frame — index 0 along
:math:`\vec Q`, index 1 in-plane perpendicular, index 2 vertical, index 3
energy. Projecting it into an ``hkle`` frame is a separate step handled by
:class:`tavi.library.resolution.ellipsoid.ResolutionEllipsoid`:

.. math::

   M_{hkle} = C^T M_{local} C, \qquad
   C_{3\times3} = 2\pi \, R_{lab \leftarrow local} \, R \, UB \, W

where :math:`R` is the sample rotation matrix from the goniometer angles and
:math:`W = [w_1\ w_2\ w_3]` holds the user's projection vectors.

Widths along a chosen axis :math:`i` come in two flavours:

.. math::

   \mathrm{FWHM}_{coh} = \frac{2\sqrt{2\ln 2}}{\sqrt{M_{ii}}}

is a *cut* through the ellipsoid (``coh_fwhm``), while the incoherent width
(``incoh_fwhm``) integrates the other three axes out with a quadric projection
first and is therefore broader.

The full argument-level description — frames, axis ordering, and how the
rotation matrix is obtained — is in :ref:`resolution`.

References
==========

The papers the implementation follows are collected in ``docs/papers/``:
Cooper & Nathans (1967), Werner & Pynn (1971), Stoica & Popovici (1975),
and Mitchell (1984).
