.. _theory:

Theoretical considerations
##########################

Neutron polarization techniques are critical tools for separating the spin-dependent and nuclear scattering components in a diverse range of quantum materials. In the current suite of direct geometry spectrometers (DGS) at the facility, HYSPEC is the only instrument specifically equipped to perform polarized neutron experiments on both powder and single crystal samples.
The underlying physics of polarized neutron scattering stems from the different interactions between neutrons and the sample, depending on the relative orientation of the neutron spin polarization and the momentum transfer vector, :math:`\vec Q`. Specifically, certain types of magnetic interactions become visible when the neutron polarization :math:`\vec P` is aligned parallel to the momentum transfer, while other types of interactions are highlighted when these two directions are perpendicular.

In the time of flight instruments, the direction of momentum transfer varies not only from detector to detector, but also as a function of the incident and scattered neutron energies. Therefore the measurements contain points witm many different angles between :math:`\vec P` and  :math:`\vec Q`.

The planning tool assumes that the polarization angle of the neutron with respect to the laboratory frame is given by the user. One can then use the kinematic/ dynamic equations of neutron scattering to determine the angle between the momentum transfer :math:`\vec Q` and the laboratory frame. In this section, we are going to use the Mantid reference frame convention, with the incident beam along :math:`\hat z`, vertical direction :math:`\hat y`, and :math:`\hat x` in the horizontal plane. Also, since HYSPEC is a vertically focusing instrument, we can ignore the :math:`y` components.

.. math::

    \begin{align}
        E_{i,f}&=\frac{\hbar^2 k_{i,f}^2}{2m_n}\\
        \Delta E&= E_i-E_f\\
        \vec Q&=\begin{pmatrix}-k_f\sin\theta\\k_i-k_f\cos\theta
        \end{pmatrix}\\
        Q^2&=k_i^2+k_f^2-2k_ik_f\cos\theta
    \end{align}

Here :math:`E_i` and :math:`E_f` are incident and final energies, :math:`k_{i,f}` are the initial and final momentum, :math:`\Delta E` is the energy transfer, :math:`\vec Q` is the momentum transfer, and :math:`\theta` is the scattering angle, which we assume that is the in-plane component only.

From here,

.. math::

    \cos\theta=\frac{k_i^2+k_f^2-Q^2}{2k_i k_f}

Note that scattering can happen both left and right of the incident beam direction. Both have the same magnitude of the scattering angle, but the sign is different. Since this would make difficult the determination of the Scharpf angle, the experimental workflow would place detectors on only one side. This allows us to determine the sign of the scattering angle.

Denoting the polarization direction :math:`\vec P=\begin{pmatrix}\sin\alpha_P\\\cos\alpha_P\end{pmatrix}`, one can get the Scharpf angle direction using

.. math::

    \cos\alpha_s=\frac{\vec P\cdot \vec Q}{|Q|}

For both powder and single crystal experiments, the user will want to find :math:`\alpha_s` at a given energy and momentum transfer. Since different directions of the same momentum transfer in the crystal are achieved by rotating the sample around the laboratory frame, it is enough to plot the Scharpf angle in a powder mode, and overplot the magnitude of :math:`\vec Q` given the lattice parameters.

.. math::

    |Q|=2\pi\left|B\begin{pmatrix}
        h\\k\\l
    \end{pmatrix}\right|

The form of the :math:`B` matrix will be the one used in Mantid (https://docs.mantidproject.org/nightly/concepts/Lattice.html).

For simplicity, sometimes it's useful to see the angle between the :math:`\vec Q` and the beam direction. Using the cosine law,

.. math::

    \cos\angle(\vec Q,\hat z)=\frac{k_i^2+Q^2-k_f^2}{2k_i Q}

Since we can measure on both sides of the incident beam, the sign of this angle is taken to be opposite of the sign of the detector tank angle.
