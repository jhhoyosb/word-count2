.. raw:: latex

   \maketitle 

Introduction
============

The study of plasmas by analytical methods alone is not
sufficient.:raw-latex:`\cite{BirdsallLangdon,tutorial}` In many
situations plasmas can be considered to be collisionless, particularly
in space and astrophysical environments, where short-range binary
interactions can be neglected.:raw-latex:`\cite{Treumann}` An approach
that uses the collisionless approximation is the particle-in-cell (PIC)
method, which was proposed in the 1950s and has become widely used due
to the great advances in computing power. The method was developed by
John Dawson and Oscar Buneman and was first applied to systems of
:math:`100` to :math:`1000` particles.:raw-latex:`\cite{Tskhakaya2008}`
The method was extended by Dawson:raw-latex:`\cite{Dawson}` in 1983. In
1988 Hockney and Eastwood published a highly acclaimed
book,:raw-latex:`\cite{Hockney_1988}` and in 1991 Birdsall and Langdon
published the most widely recognized book on the
method.:raw-latex:`\cite{BirdsallLangdon}` Both books describe the PIC
method and discuss different simulations and analyses.

Simulations are widely used to study different kinds of
microinstabilities and waves common in space plasmas, such as the
electron-electron two-stream instability,:raw-latex:`\cite{Hutch}`
Buneman instability,:raw-latex:`\cite{Rajawat_2017}` and ion-acoustic
waves.:raw-latex:`\cite{Koen}` There are now open source programs such
as Jupyter notebooks for educational plasma
simulations,:raw-latex:`\cite{Ucla}` and the graphical user interface
developed by Omura and Matsumoto.:raw-latex:`\cite{Omura}`

Reference  discusses the implementation of the two-dimensional
electrostatic PIC method in detail. Although two dimensions yields a
more realistic view of plasmas, the computational cost and the
complexity of the simulations is much greater than one dimension. As
will be seen in Sec. II, PIC codes are usually built over an Eulerian
mesh, and two dimensions usually requires the order of
:math:`32 \times 32 = 1024` grid points, in contrast to one dimension
(1d) for which the minimum is the order of :math:`32` grid points,
meaning that the computational cost increases dramatically in 2d
systems.:raw-latex:`\cite{BirdsallLangdon}` Moreover, collisionless
plasmas contain many particles in a characteristic volume known as the
Debye sphere, whose radius is the Debye length :math:`\lambda_D`. The
number of particles in this volume is :math:`N_D\sim \lambda_D^d`. If a
3d system has :math:`N_D= 10^6` for example, then
:math:`N_D\approx (10^6)^{2/3}=10^4` in 2d, and
:math:`N_D\approx (10^6)^{1/3}=10^2` in 1d, which implies a lower
computational effort for the 1d
system.:raw-latex:`\cite{BirdsallLangdon}` Understanding two-dimensional
simulations requires a more advanced mastery of the PIC method and of
the physics involved, and thus a one-dimensional system provides a good
starting point for students.

Although Refs.  explain the PIC method in detail, its understanding
requires the mastery of numerical methods and advanced electromagnetic
concepts. The main objective of this paper is to explain the PIC method
so that students who have taken a course in integral calculus and a
basic course in electrostatics can understand the method. In this
spirit, we consider the simplest plasma, i.e., a Coulomb plasma, for
which the system is free of magnetic fields and the plasma waves contain
only electric fields varying in a longitudinal direction. A graphical
user interface (PlasmAPP) is developed to help teachers and students
study basic plasma phenomena in an interactive way.

Our numerical implementation of the PIC model is based on the Lapenta
open source code,:raw-latex:`\cite{LapentaGIT}` for which the plasma is
composed of moving electrons and ions fixed in the background. We
extended the code to include the motion of ions, two interpolation
methods, different numerical methods for the solution of the
differential equations, and multiple diagnostics. The numerical methods
explained in Sec. III are the direct integration of Gauss’s law of
electrostatics using the trapezoid method to approximate integrals,
Euler’s method for the time integration of Newton’s second law for
particle dynamics, and the nearest grid point as the simplest
interpolation method. More advanced numerical methods, which can
increase the efficiency and the speed of the simulations, are explained
in the supplemental
material.:raw-latex:`\cite{PlasmAPPGIT, PlasmAPPAJP}`

To illustrate the application of the one-dimensional PIC method, we
discuss simulations of the electron-electron two-stream instability, a
well-known phenomena corresponding to an initial condition where two
groups of electrons move in opposite directions. This condition departs
from thermal equilibrium characterized by a single Maxwellian
distribution function.:raw-latex:`\cite{REI65}` Reference  also studied
this phenomenon using a two-dimensional simulation. They simulated two
electron beams moving with counter-propagating initial velocities along
a single Cartesian direction. The subsequent evolution of the system in
search of thermal equilibrium leads to the appearance of electric fields
pointing mainly along this Cartesian direction. However, it is
sufficient to study this situation with a one-dimensional model as we do
in the present work. Our code includes comments to facilitate
understanding and is freely
available.:raw-latex:`\cite{PlasmAPPGIT,PlasmAPPAJP}` Interested readers
are guided to choose parameters that can be used to simulate other
plasma phenomena such as the Buneman instability, ion-acoustic waves,
and Langmuir waves.

The paper is organized as follows. Section II explains the basic
ingredients of one-dimensional electrostatic PIC simulations.
Section III discusses the PIC method and explains its numerical
implementation. Section IV discusses the structure of the graphical user
interface. Section V contains the results and analysis of the
electron-electron two-stream instability simulations. Section VI
discusses our conclusions and Sec. VII includes some suggested problems
for interested readers.

Theoretical background
======================

Because plasmas are composed of millions of particles, the dynamics of
each charge as a function of the self-generated electric forces must be
followed to find the temporal evolution of the system. This approach
requires finding the force experienced by each particle produced by the
remaining particles and solving Newton’s second law for all of the
charges. For a plasma composed of :math:`N` particles, :math:`N-1`
forces must be calculated to obtain the force experienced by a single
charge. If this process is repeated for each of the :math:`N` particles,
:math:`N(N-1)/2` forces must be calculated, which requires a high
computational cost.:raw-latex:`\cite{BirdsallLangdon,tutorial,Lapenta3}`
It is therefore necessary to use methods that approximate the behavior
of the plasma at a lower computational cost. There are two general ways
to achieve this task. One is the plasma fluid approach, which studies
the evolution of macroscopic or bulk properties such as the density,
pressure, mean velocities, and mean energy.:raw-latex:`\cite{Ledvina}`
This approach is useful for understanding plasma dynamics for phenomena
occurring on large spatial and temporal scales, in contrast with the
small scales associated with the discrete nature of
particles.:raw-latex:`\cite{goldston}` Microscopic temporal and spatial
scales are best modeled by a kinetic approach that considers the
discrete nature of particles in constant motion due to the
self-consistent electromagnetic fields. Because a detailed kinetic
description is not computationally efficient due to the large number of
particles, it is necessary to use a reduced but sufficiently detailed
model that contains information about the particle velocity distribution
function.:raw-latex:`\cite{Hockney_1988}`

A kinetic description of plasmas can be achieved by using the Vlasov
equation, which gives the spatio-temporal evolution of the distribution
function of a single particle.:raw-latex:`\cite{Bittencourt}` This
function is defined in phase space consisting of the spatial coordinates
and particle velocities of the particles.:raw-latex:`\cite{bellan_2006}`
The distribution function gives the probability of finding a particle in
a small region of phase space.:raw-latex:`\cite{REI65}` The Vlasov model
corresponds to a mean-field theory that evolves the particles by
long-range electromagnetic fields calculated as a statistical average
over the ensemble of particles. This model is a good approximation for
time scales smaller than the characteristic time of the collisions and
for large length scales compared to the average particle
separation.:raw-latex:`\cite{krall1973principles}` If the
self-consistent coupling of the Vlasov equation with Maxwell’s equations
is considered, the set of Vlasov-Maxwell equations is
obtained.:raw-latex:`\cite{Chaudhary18}`

An approximate solution of the Vlasov equation can be found by finite
differences on an Eulerian mesh. During the evolution, particles occupy
only a finite region of phase space.:raw-latex:`\cite{RUHL1995388}` For
this reason, the implementation of Vlasov codes requires the
implementation of arrays of which only a small part is used at any one
time, which implies an inefficient use of memory and
processing.:raw-latex:`\cite{pukhov}`

An approach with low computational cost is the PIC method. The system
consists of superparticles, which represent the dynamics of many, but
not too many, real particles moving in continuous phase space
(Lagrangian description). The evolution of the distribution function of
each species is obtained from the positions and velocities of these
superparticles.

If we fix one particle in the plasma, the other particles will move to
reduce the electric potential produced by the fixed particle. The
resulting shielded electric potential will decay exponentially with a
length scale called the Debye length,:raw-latex:`\cite{Bittencourt}`
which is defined as

.. math::

   \lambda_{D,\,\alpha}=\frac{v_{\rm th,\,\alpha}}{\omega_{p,\,\alpha}
   },

 where
:math:`\omega_{p,\,\alpha}^2= n_\alpha q_\alpha^2/m_\alpha\epsilon_0` is
the plasma frequency,
:math:`v_{\rm th,\,\alpha}^2=k_B T_\alpha/m_\alpha` is the thermal
velocity, :math:`n_\alpha` is the number density, :math:`q_\alpha` the
charge, :math:`m_\alpha` the mass, and :math:`T_\alpha` the temperature
of the particle velocity distribution function for species
:math:`\alpha`. Although the electric potential produced by one particle
decays with distance, particles closer to the edge of the Debye sphere
of the fixed particle can interact with particles nearby but outside the
Debye sphere, producing a collective interaction that can propagate
through the plasma. Hence, superparticles can be constructed as finite
clouds of particles subject to electrostatic interactions with other
superparticles that are far away, but experience a weak interaction when
they are close to each other.:raw-latex:`\cite{Dawson}` The behavior of
the spatial density is not modeled in regions smaller than the cloud
size. In this way it makes sense to construct the charge density and
solve for the electric force from the self-consistent electric field
through Maxwell equations which are solved on a fixed Eulerian mesh in
configuration space.:raw-latex:`\cite{pukhov}` The field can then be
interpolated at the position of each superparticle to calculate their
acceleration.:raw-latex:`\cite{Navarro}`

Electrostatic one-dimensional PIC algorithm
===========================================

In this work we consider the simplest configuration which is a
one-dimensional system. Many of the technical subtleties of :math:`2`\ d
simulations are already present in :math:`1`\ d simulations, both in
physics and implementation, so before tackling the task of constructing
a simulation for :math:`2`\ d plasmas, it is convenient to understand
the difficulties that already appear in :math:`1`\ d plasma simulations.
In this type of simulation, harmonic waves with wave propagation vector
:math:`\vec{k}` and frequency :math:`\omega` have their electric field
parallel to :math:`\vec{k}` (parallel propagation). Therefore, if one
then wants to advance to electromagnetic or magnetized problems, it is
useful to start with :math:`1`\ d and then advance to :math:`2`\ d
(oblique propagation).:raw-latex:`\cite{BirdsallLangdon}` There are a
number of current problems in basic and applied research that can be
studied with :math:`1`\ d PIC
simulations.:raw-latex:`\cite{Sania,Pandey}`

As mentioned, the PIC method is based on the use of superparticles whose
positions and velocities follow a continuous path in phase space, while
spatial macroscopic quantities, such as the charge density and
self-generated electric field, are calculated at discrete spatial points
of a grid (Eulerian description). Figure \ `[fig:mesh] <#fig:mesh>`__
shows how the Eulerian space grid is created for one dimension. The
solid line represents the length of the system ranging from :math:`0` to
:math:`L`. The grid is divided into :math:`N_g` nodes denoted by
discrete indices :math:`g = 1,\; 2, \;3,\;\ldots,\;N_g`. The system is
periodic, so that the simulation box of length :math:`L` represents a
section of an infinite plasma that repeats periodically in space. For a
periodic system, nodes divide the total length :math:`L` into
:math:`N_g` cells of size :math:`\Delta x`, centered on the grid points
as shown in Fig. \ `[fig:mesh] <#fig:mesh>`__. The grid node positions
are :math:`x_g=g\Delta x`, and the length of the system is
:math:`L=N_g\Delta x`.:raw-latex:`\cite{Hockney_1988}`

.. raw:: latex

   \centering

Periodic boundary conditions are used so that a spatially variable
quantity :math:`G(x)` has a spatial period :math:`L`, i.e.,
:math:`G(x+L)=G(X)`. Thus, the last point of the mesh :math:`N_g` is the
same as the first point :math:`g=0` (which in not part of the array,
which starts at :math:`g=1`), and the point :math:`N_g+1` is the same
point as :math:`g=1`. If a superparticle moves outside the simulation
box, it returns to a position inside the
box.:raw-latex:`\cite{BirdsallLangdon}` The distance between nodes
should be smaller than the Debye length, and within each cell there are
a large number of superparticles so that we can observe collective
effects.

The PIC algorithm is depicted in
Fig. \ `[fig:PIC_ALGO] <#fig:PIC_ALGO>`__ and is summarized in the
following.

.. raw:: latex

   \centering

.. raw:: latex

   \tikzset {_03gz55301/.code = {\pgfsetadditionalshadetransform{ \pgftransformshift{\pgfpoint{89.1 bp } { -128.7 bp }  }  \pgftransformscale{1.32 }  }}}

.. raw:: latex

   \pgfdeclareradialshading{_4se5zmyh3}{\pgfpoint{-72bp}{104bp}}{rgb(0bp)=(0.96,0.65,0.14);
   rgb(0bp)=(0.96,0.65,0.14);
   rgb(25bp)=(0.82,0.01,0.11);
   rgb(400bp)=(0.82,0.01,0.11)}

.. raw:: latex

   \tikzset {_8bltk9cid/.code = {\pgfsetadditionalshadetransform{ \pgftransformshift{\pgfpoint{89.1 bp } { -128.7 bp }  }  \pgftransformscale{1.32 }  }}}

.. raw:: latex

   \pgfdeclareradialshading{_pfrxcn3mf}{\pgfpoint{-72bp}{104bp}}{rgb(0bp)=(0.96,0.65,0.14);
   rgb(0bp)=(0.96,0.65,0.14);
   rgb(25bp)=(0.82,0.01,0.11);
   rgb(400bp)=(0.82,0.01,0.11)}

.. raw:: latex

   \tikzset {_epa2n6nio/.code = {\pgfsetadditionalshadetransform{ \pgftransformshift{\pgfpoint{89.1 bp } { -128.7 bp }  }  \pgftransformscale{1.32 }  }}}

.. raw:: latex

   \pgfdeclareradialshading{_akc42el55}{\pgfpoint{-72bp}{104bp}}{rgb(0bp)=(0.96,0.65,0.14);
   rgb(0bp)=(0.96,0.65,0.14);
   rgb(25bp)=(0.82,0.01,0.11);
   rgb(400bp)=(0.82,0.01,0.11)}

.. raw:: latex

   \tikzset {_97etwtbmz/.code = {\pgfsetadditionalshadetransform{ \pgftransformshift{\pgfpoint{89.1 bp } { -128.7 bp }  }  \pgftransformscale{1.32 }  }}}

.. raw:: latex

   \pgfdeclareradialshading{_q9hxssdjf}{\pgfpoint{-72bp}{104bp}}{rgb(0bp)=(0.96,0.65,0.14);
   rgb(0bp)=(0.96,0.65,0.14);
   rgb(25bp)=(0.82,0.01,0.11);
   rgb(400bp)=(0.82,0.01,0.11)}

[fig:PIC_ALGO]

#. Initialize the distribution function of the superparticles indexed by
   :math:`p` belonging to species
   :math:`\alpha = \text{ions or electrons}`, i.e., assign initial
   positions and velocities to the superparticles. For simplicity, the
   superparticles are initially positioned uniformly in space and if
   desired, a small sinusoidal perturbation is added. The initial
   velocities are assigned according to a Maxwellian distribution
   corresponding to a system in thermal
   equilibrium.:raw-latex:`\cite{REI65}`

   The width of the Maxwellian is determined by
   :math:`v_{\rm th,\,\alpha}`. The function is centered on the speed
   :math:`v_0`, which is determined by the nature of the problem being
   studied.

#. Calculate the charge density at the grid points :math:`\rho_g`.
   Because the coordinate :math:`x_{p,\,\alpha}` of the superparticles
   takes on continuous values, it is necessary to determine the
   contribution of each superparticle to the charge density at the grid
   points. To do so, we define an interpolation function which
   determines the weight that each superparticle contributes to the
   physical properties at the grid points. The simplest method is called
   the nearest grid point. If the distance between the grid point and
   the center of superparticle :math:`|x_g-x_{p,\,\alpha}|` is less than
   :math:`\Delta x/2`, a weight of :math:`1` is assigned to that grid
   point. This method can be expressed using the b-spline flat-top
   function :math:`b_0` as

   .. math::

      b_0(\chi) =
      \begin{cases}
          1 & |\chi|<1/2 \\
          0 &  \mbox{otherwise}.
          \label{b-spline}
          \end{cases}

    The corresponding interpolation function for species :math:`\alpha`
   is given by

   .. math::

      W_\alpha(x_g-{x}_{p,\,\alpha})= b_0\left(\frac{x_g-{x}_{p,\,\alpha}}{\Delta x}\right).
          \label{interp1}

   This interpolation scheme is shown in
   Fig. \ `[fig:WNGP] <#fig:WNGP>`__. In more advanced methods the
   superparticle is assigned to more than one grid
   point.:raw-latex:`\cite{BirdsallLangdon,Hockney_1988,Supplementary}`

   .. raw:: latex

      \centering

#. The charge density at the grid points :math:`g` for species
   :math:`\alpha` is given by:raw-latex:`\cite{Lapenta3}`

   .. math::

      \rho_{g,\,\alpha} =\frac{q_{p,\,\alpha}}{\Delta x} \sum_p  W_\alpha(x_g-{x}_{p,\,\alpha}).
          \label{e1}

    The charge of the superparticle is defined as
   :math:`q_{p,\,\alpha}=q_\alpha N_{r,\,\alpha}` where
   :math:`N_{r,\,\alpha}` is the number of real particles in a
   superparticle.

#. The total charge density at grid point :math:`g` is

   .. math::

      \rho_g = \sum_\alpha \rho_{g,\,\alpha}.
             \label{rhogrid}

#. Given the charge density at the grid points, we calculate the
   electric field from Poisson’s equation at these nodes.

   .. math::

      \frac{d^2\phi(x)}{dx^2}=-\frac{\rho(x)} {\epsilon_0}.
              \label{poisson-}

   If we integrate both sides of Eq. \ `[poisson-] <#poisson->`__ from
   the left boundary :math:`x=0` to an arbitrary point :math:`x`, we
   obtain

   .. math::

      E(x) = E(0) + \frac{1}{\epsilon_0}\int_0^x \rho(x')dx',
              \label{continuousE}

    where we have used the relation :math:`E(x)=-d\phi(x)/dx`. Note that
   Eq. \ `[continuousE] <#continuousE>`__ requires the value of the
   electric field at the left boundary :math:`E(0)`. Because of periodic
   boundary conditions, we also have :math:`E(L)=E(0)`. To obtain an
   expression for :math:`E(0)` in terms of the charge density, we
   integrate Eq. \ `[poisson-] <#poisson->`__ twice and use
   :math:`\phi(0)=\phi(L)` to obtain

   .. math::

      E(0)=-\frac{1}{L}\int_0^L \left(\int_0^{x'} \rho(x'') dx'' \right) dx',
           \label{E0}

    Because the charge density is evaluated at the grid points, it is
   necessary to compute the integrals in
   Eqs. \ `[continuousE] <#continuousE>`__ and `[E0] <#E0>`__ using the
   values of the electric field. For this purpose, the trapezoidal rule
   is used, so that the electric field at each grid point is given
   by:raw-latex:`\cite{Burden1989}`

   .. math::

      E_g = E_{g=0} + \frac{\Delta x}{\epsilon_0}\sum_{j=1}^{g-1}\left(\frac{\rho_{j}+\rho_{j-1}}{2}\right),
              \label{Direct_Eg}

    where :math:`\rho_0=\rho_{N_g}`. In the same way, it is possible to
   obtain the electric potential from the relation

   .. math::

      \phi(x) = \phi(0) - \int_0^x E(y) dy.
              \label{phix}

    Because the reference potential :math:`\phi(0)` in
   Eq. \ `[phix] <#phix>`__ can be selected arbitrarily, we set
   :math:`\phi(0)=0` for simplicity. We evaluate the integral in
   Eq. \ `[phix] <#phix>`__ numerically and express :math:`\phi_g` as

   .. math:: \phi_g =  \phi_{g=0} - \Delta x \sum_{j=1}^{g-1}\left(\frac{E_{j}+E_{j-1}}{2}\right),

    where :math:`\phi_{g=0}=0=\phi_{N_g}`.

#. To find the electric field :math:`E_{p,\,\alpha}` experienced by a
   superparticle at the continuous position :math:`x_{p,\,\alpha}`, the
   nearest grid point method (see Fig. \ `[fig:NGP_E] <#fig:NGP_E>`__)
   is used so that a superparticle in a cell centered at :math:`g`
   (closest grid point to the superparticle) interacts with the electric
   field :math:`E_g`. That is, the same interpolation function used in
   Eq. \ `[interp1] <#interp1>`__ is applied. Therefore,
   :math:`E_{p,\,\alpha}` is given by

   .. math::

      E_{p,\,\alpha}= \sum_g E_g W_\alpha(x_g-x_{p,\,\alpha}).
              \label{Epalpha}

    The same expression can be applied to more general interpolation
   methods as shown in Refs. . The electric force on the superparticle
   is

   .. math::

      F_{p,\,\alpha} = q_{p,\,\alpha}E_{p,\,\alpha}.
          \label{Lorentz}

   .. raw:: latex

      \centering

   [fig:NGP_E]

#. The dynamics of the center of mass of the superparticles is governed
   by Newton’s equations:raw-latex:`\cite{Hockney_1988,Lapenta3}`

   .. math::

      \begin{aligned}
          \frac{d{x}_{p,\,\alpha}}{dt} & ={v}_{p,\,\alpha},
          \label{posicionC} \\
              \frac{d{v}_{p,\,\alpha}}{dt} & =\frac{F_{p,\,\alpha}}{m_{p,\,\alpha}}.
              \label{velocidadC}
          \end{aligned}

    The solution of Eqs. \ `[posicionC] <#posicionC>`__ and
   `[velocidadC] <#velocidadC>`__ are obtained numerically using the
   Euler algorithm

   .. math::

      \begin{aligned}
          {x}_{n+1,\,p,\alpha} & = {x}_{n,\,p,\alpha} +{v}_{n,\,p,\alpha}\Delta t,
          \label{Euler_posición} \\
          {v}_{n+1,\,p\alpha} & ={v}_{n,\,p\alpha}+ \frac{q_{p,\,\alpha}}{m_{p,\,\alpha}}E_{n,\,p,\alpha} \Delta t.
          \label{Euler_velocidad}\end{aligned}

    where :math:`t_n=n \Delta t` and :math:`\Delta t` is the time step.

As can be seen in Eq. \ `[Euler_velocidad] <#Euler_velocidad>`__, the
acceleration depends on the charge-to-mass ratio of the superparticles.
Because the charge and mass of a superparticle :math:`p` of species
:math:`\alpha` is given by :math:`q_{p,\,\alpha}=N_{r,\,\alpha}q_\alpha`
and :math:`m_{p,\,\alpha}=N_{r,\,\alpha}m_\alpha`, where
:math:`N_{r,\,\alpha}` is the number of real particles of species
:math:`\alpha`, we have:raw-latex:`\cite{pukhov}`

.. math:: \frac{q_\alpha}{m_\alpha}=\frac{q_{p,\,\alpha}}{m_{p,\,\alpha}}.

 Hence, the center of mass of a superparticle follows the same
trajectory as a real particle because they experience the same
acceleration.

PlasmAPP
========

To make plasma simulations more accessible, we created a program with a
graphical user interface which we call PlasmAPP. The program allows
users to choose the Euler, leapfrog, or fourth-order Runge-Kutta
algorithm for the solution of the equations of
motion.:raw-latex:`\cite{Burden1989,BirdsallLangdon}` Users can explore
the advantages and disadvantages of these algorithms in terms of
numerical stability, speed, and numerical accuracy. The method discussed
in Sec. III to find the electric field was the direct integration of
Gauss’s law. The more advanced possibilities implemented in PlasmAPP
include the finite difference method and the fast Fourier
transform.:raw-latex:`\cite{Burden1989,Supplementary}` Users can also
select the interpolation method used to assign the charge to the grid
and the force calculations. Two possibilities are the nearest grid point
explained in Sec. III and the cloud-in-cell
method.:raw-latex:`\cite{BirdsallLangdon,Hockney_1988,Supplementary}`
PlasmAPP also allows the simulation of phenomena involving ion mobility
as well as a fixed neutralizing background. The diagnostics can be
chosen and the results can be saved for post-processing. The code can be
found at Refs.  and additional information of the code can be found at
Ref. .

Results
=======

We discuss the results of the electron-electron two-stream instability
using the Euler method for the equations of motion and direct
integration for the Gauss’s law. The interpolation method used is the
nearest grid point. These methods have a greater numerical error than
more advanced methods, but they were discussed here because they are
easier for students to understand who only have a knowledge of integral
calculus and basic electrostatics.

It is convenient to normalize the variables in terms of characteristic
times and lengths which define the corresponding scales of interest and
avoid the use of very small or very large numbers that lead to numerical
errors. By means of an inverse procedure, it is possible to find the
physical quantities and their real values. We measured time in terms of
the electron plasma period (:math:`\omega_{p,e}^{-1}`), lengths in terms
of the Debye length :math:`\lambda_{D,e}`, the mass in terms of the
electron mass :math:`m_e`, and the charge in terms of the proton charge
:math:`e`. The vacuum permittivity :math:`\epsilon_0` and the Boltzmann
constant :math:`k_B` are set to one. The other units of the physical
quantities such as the velocity, electric potential, electric field, and
energy among others, are derived from these fundamental units.

Electron-electron two-stream instability
----------------------------------------

The electron-electron two-stream instability consists of two electron
beams with opposite velocities in an immobile ion background. The
velocity distribution function is given
by:raw-latex:`\cite{Bland_n_2017}`

.. math::

   f({v}_x)  = \frac{n_0}{2\sqrt{2\pi}v_{\rm th,e}}e^{\frac{-(v-v_{0})^2}{2v_{\rm th,e}^2}} \nonumber  + \frac{n_0}{2\sqrt{2\pi}v_{\rm  th,e}}e^{\frac{-(v+v_{0})^2}{2v_{\rm th,e}^2}}.
       \label{dist2}

To simulate the electron-electron two-stream instability, we used the
parameters :math:`L=64`, number of time steps :math:`N_t = 16000`,
:math:`\Delta t= 0.1`, number of grid points :math:`N_g = 256`, number
of electrons of the two beams :math:`N_{e,1}=N_{e,2}=10000`, number of
ions fixed to the background :math:`N_{f}=20000`, :math:`v0_{e,1}=5`,
:math:`v0_{e,2}=-5`, thermal speed of the two beams
:math:`v_{\rm th,e,1}=v_{\rm th,\,e,2} = 1`, charge-to-mass relation of
the beams :math:`r_{e,1}=r_{e,2}=-1`, and the plasma frequency
:math:`\omega_{p,e}=1`.

The initial state of the system is shown in
Fig. \ `[fig:T_PS_DF_it1] <#fig:T_PS_DF_it1>`__, where the phase space
and the distribution function is displayed. The latter shows the
presence of two peaks centered at :math:`v0_{e,1}` and :math:`v0_{e,2}`,
which correspond to five times the thermal speed. These peaks represent
the two electron beams. The initial state is charge neutral, and
therefore the electric field in the grid is zero everywhere.

.. raw:: latex

   \centering

.. raw:: latex

   \centering

.. figure:: 5a.eps
   :alt: 
   :name: fig:TPS1

.. raw:: latex

   \centering

.. figure:: 5b.eps
   :alt: 
   :name: fig:TDF1

[fig:T_PS_DF_it1]

The initial particle velocity distribution function eventually generates
an electric field, allowing the beams to interact. Initially the
electric field perturbations are waves that grow exponentially in time,
signaling an instability (see
Fig. \ `[fig:T_PS_Phi_E_it105] <#fig:T_PS_Phi_E_it105>`__).

.. raw:: latex

   \centering

.. raw:: latex

   \centering

.. figure:: 6a.eps
   :alt: 

.. raw:: latex

   \centering

.. figure:: 6b.eps
   :alt: 

.. raw:: latex

   \centering

.. figure:: 6c.eps
   :alt: 

[fig:T_PS_Phi_E_it105]

.. raw:: latex

   \centering

.. raw:: latex

   \centering

.. figure:: 7a.eps
   :alt: 

[T_PS_3500]

.. raw:: latex

   \centering

.. figure:: 7b.eps
   :alt: 

[T_Phi_3500]

.. raw:: latex

   \centering

.. figure:: 7c.eps
   :alt: 

[fig:T_Field_3500]

[fig:T_PS_Phi_E_it3500]

When the particles interact, they start to form electron phase space
holes as observed in
Fig. \ `[fig:T_PS_Phi_E_it3500] <#fig:T_PS_Phi_E_it3500>`__, which
corresponds to the system at time :math:`t=350`. As can be seen in
Fig. \ `[T_Phi_3500] <#T_Phi_3500>`__, the potential reaches a maximum
at the center of the electron hole. Because the electric field is equal
to the negative derivative of the electric potential, this maximum
separates the two regions where the electric field changes sign (see
Fig. \ `[fig:T_Field_3500] <#fig:T_Field_3500>`__). Because the
direction of the electric force is opposite to the direction of the
electric field due to the sign of the charge, particles that are to the
right (left) of the maximum will experience forces to the left (right),
so particles will tend to change their direction of motion, which
implies a curvature in phase space that is visualized as a hole (see
Fig. \ `[T_PS_3500] <#T_PS_3500>`__).:raw-latex:`\cite{Hutchinson}`
Similar results were obtained in two dimensions for simulations of the
electron two-stream instability with fixed ions in the presence of a
constant external magnetic field.:raw-latex:`\cite{Jaime}` In their
simulation they used the same initial condition as
Eq. \ `[dist2] <#dist2>`__ with the two electron beams moving in
opposite directions along the :math:`x`-direction (see Fig. 5a of
Ref. ). They showed that the initial value of the electric potential
:math:`\phi` is zero, and as the system evolves, waves that vary mainly
along the :math:`x` coordinate begin to appear (see Fig. 6 of Ref. ) as
the system attempts to reach thermal equilibrium. It was also shown that
the electric potential is positive at localized regions where electron
holes in phase space are formed (see Figs. 5 and 6 of Ref. ). Because
both electron beams in Ref.  have initial velocities in just one
direction (:math:`x`), the spatial variation of :math:`\phi` through the
:math:`y` coordinate is not significant. Also, the presence of the
constant external magnetic field does not influence the evolution in the
electrostatic limit because currents are not generated along the
:math:`x` direction according to Ampere’s law. :raw-latex:`\cite{Jaime}`
Thus, our 1d simulation contains the relevant physics of the instability
for this particular initial condition, with the advantage of less
computational effort in comparison to the 2d simulations. However, the
2d code implemented by Ref.  can be used to simulate the more general
oblique propagation.

Figure \ `[fig:T_Energy_152] <#fig:T_Energy_152>`__ shows how the energy
behaves in the time interval :math:`t=0` to :math:`t=20`. We observe
that the kinetic energy :math:`\varepsilon_k` (dashed-dotted line) is
much larger than the potential energy :math:`\varepsilon_p` (solid
line). When :math:`\varepsilon_k \gg\varepsilon_p`, the system
corresponds to a weakly coupled plasma, where long-range interactions
dominate over short-range collisions.

The movement of the electron beams with respect to each other causes
density perturbations due to the electric force between the particles.
These perturbations cause a linear instability characterized by an
initial temporal exponential growth of the electric field according to
the linear instability theory in
plasmas.:raw-latex:`\cite{Chen,Aper,BirdsallLangdon}` This behavior is
shown in Fig. \ `[fig:T_Energy_152] <#fig:T_Energy_152>`__ where the
initial growth of the potential energy occurs for
:math:`10\lesssim t \lesssim 14`. At :math:`t=15.2`, the instability is
first saturated and the dynamics becomes nonlinear as shown in
Fig. \ `[fig:T_PS_152] <#fig:T_PS_152>`__. After a long time, the beam
interaction eventually causes the particles’ speeds to converge to the
Maxwell-Boltzmann distribution. This long time convergence can be
observed using the more advanced methods discussed in Refs. .

.. raw:: latex

   \centering

.. raw:: latex

   \centering

.. figure:: 8a.eps
   :alt: 

[fig:T_Energy_152]

.. raw:: latex

   \centering

.. figure:: 8b.eps
   :alt: 

[fig:T_PS_152]

[fig:T_PS_Ener_it152]

The relative percentage error of the total energy at each iteration with
respect to the initial total energy was calculated and found to be
:math:`0.57\%`, indicating that the energy is conserved with reasonable
accuracy. An analogous error analysis was performed for the total
momentum, which should be zero because the two electron beams have the
same density and the same speed but move in opposite directions. The
maximum percentage error of the momentum is
:math:`5.54\times10^{-10}\,\%`.

Discussion
==========

To understand the behavior of systems as complex as plasmas, it is
helpful to use the simplest possible tools. Visualizing the complex
behavior of plasmas provides additional insights, particularly behavior
that is not easily accessible by theory.

Plasma evolution involves self-consistent dynamics between the
electromagnetic fields and particles, because particles are the source
of these fields and particle trajectories are determined by the fields.
We have discussed the particle-in-cell method using direct and simple
numerical methods and focused on the implementation of a simplified
plasma system in the electrostatic regime with a one-dimensional
geometry where all physical quantities vary along a single direction.

As an example of the particle-in-cell method and the PlasmAPP program,
we simulated the instability of two electron streams where the formation
of electron holes in phase space is observed, leading to the generation
of bipolar pulses of the electric field, which is consistent with
measurements in astrophysical environments.:raw-latex:`\cite{Picket}`
Similar behavior is obtained in one and two dimensions if the initial
condition corresponds to electron beams moving along a single Cartesian
coordinate. This example shows that the behavior of a system can be
simulated by a simpler one-dimensional simulation without the need to
introduce the complexities of a simulation in higher dimensions.

Readers can download PlasmAPP from Refs.  and perform simulations of
other phenomena of interest, including ion motion, and explore other
numerical methods. We hope that we have provided a guide for readers
beginning the study of plasma physics and kinetic simulations.

Suggested problems
==================

*Problem 1: Write your own code*. We present a brief guide to writing
your own electrostatic one-dimensional code using the Euler method for
integrating the motion of the particles and direct integration with the
trapezoidal rule for determining the electric field. The following steps
are for an electron beam with a fixed ion background.

-  Define the parameters of the system: length of the system :math:`L`,
   the number of grid points :math:`N_g`, the cell size
   :math:`\Delta x`, number of time steps :math:`N_t`, and the time step
   :math:`\Delta t`. Then define the parameters for the superparticles:
   the number of superparticles for electrons :math:`N_{\rm e}`, the
   electron plasma frequency :math:`\omega_{p,e}`, the initial velocity
   :math:`v0_e`, the thermal velocity :math:`v_{\rm th,\,e}`, and the
   charge-to-mass ratio :math:`r_e`. The charge of the superparticle is
   given by :math:`Q_{\rm e} = \omega_{p,e}^2/[r_e (N_{\rm e}/L)]`. To
   obtain quasi-neutrality of the plasma, a background charge density
   :math:`\rho_{\rm back}` must be added
   :math:`(N_e/L)Q_{\rm e} +\rho_{\rm back} = 0`, so that
   :math:`\rho_{\rm back}= - (N_e/L)Q_{\rm e}`.

-  Place the particles equally spaced in position. One way to obtain a
   Maxwellian distribution for the velocities is to multiply the thermal
   velocity by a function that creates normally distributed random
   numbers:raw-latex:`\cite{Randn,RandnP}` and then add the initial
   velocity of the electrons.

-  Apply the interpolation function of Eq. \ `[interp1] <#interp1>`__ to
   construct the charge density at the grid points. Include periodic
   boundary conditions for :math:`x_g`, which is the position of the
   grid point closest to the superparticle. Thus, if a superparticle is
   near :math:`g=0`, the charge assignment is to the last grid point.

-  Calculate the charge density using Eq. \ `[e1] <#e1>`__, iincluding
   the background density :math:`\rho_{\rm back}`.

-  Use direct integration with the trapezoid rule to compute the
   electric field at the grid points with
   Eq. \ `[Direct_Eg] <#Direct_Eg>`__.

-  To obtain the electric field :math:`E_{p,\,\alpha}` experienced by a
   superparticle at its current position :math:`x_{p,\,\alpha}`, use
   Eq. \ `[Epalpha] <#Epalpha>`__ to interpolate this value from the
   electric fields :math:`E_g` at the grid points.

-  Calculate the new velocities of the superparticles using
   Eq. \ `[Euler_velocidad] <#Euler_velocidad>`__.

-  Use Eq. \ `[Euler_posición] <#Euler_posición>`__ to calculate the new
   particle positions. Make sure that you use periodic boundary
   conditions on the updated positions.

In the following problems use your own program or the one provided in
Refs.  to simulate the evolution of the plasma.

*Problem 2: Langmuir waves*. We can observe Langmuir waves in a plasma
in thermal equilibrium. Use the parameters :math:`L=1024`,
:math:`N_t = 8000`, :math:`\Delta t= 0.05`, :math:`N_g =8192`,
:math:`N_e=50000`, :math:`v_{\rm th,\,e}=1`, :math:`r_e=-1`,
:math:`\omega_{p,e}=1`, and set the other parameters to zero. It is
necessary in PlasmAPP to define the number of ions in the background to
be equal to the number of electrons. The dispersion relation for
Langmuir waves is:raw-latex:`\cite{HannuWeltraumBuch}`

.. math::

   \omega^2 =\pm\omega_{p,\,e}^2 \left(1+\frac{3}{2}k^2\lambda_{D,\,e}^2\right),
       \label{LangDR}

 where :math:`\omega` corresponds to the frequency associated with the
propagating harmonic fields and :math:`k` is the corresponding
wavenumber. Use your results to plot the theoretical dispersion relation
and the simulated one, and determine how they compare. To do so, apply a
fast Fourier transform in space and time to the electric field. You
should obtain a plot like the one in Fig. \ `[fig:lang] <#fig:lang>`__.

.. raw:: latex

   \centering

.. figure:: 9.eps
   :alt: Ratio of the theoretical to the simulated dispersion relation
   for Langmuir waves (see Problem 2).
   :name: fig:lang

   Ratio of the theoretical to the simulated dispersion relation for
   Langmuir waves (see Problem 2).

Vary the thermal velocity and observe how the dispersion relation
changes. If the thermal velocity becomes smaller, you will obtain cold
plasma waves with a frequency equal to :math:`\omega_{p,e}`. You will
also see that the numerical error increases because of the
non-compliance of the stability
conditions.:raw-latex:`\cite{Hockney_1988}` It is important to be aware
of the limitations of the program. Add an initial velocity to the beam
of electrons, and observe if there is a change in the dispersion
relation.

*Problem 3: Ion-acoustic and Langmuir waves*. In Problem 2 the ions were
fixed. What happens when they are mobile? Add a second species to your
program and use the parameters: :math:`L=2048`, :math:`N_t = 8000`,
:math:`\Delta t= 0.05`, :math:`N_g =8192`, :math:`N_e=8000`,
:math:`N_i=8000`, :math:`v0{_e}=0`, :math:`v0{_i}=0`,
:math:`v_{\rm th,\,e}=10`, :math:`v_{\rm th,\,i}=0`, :math:`r_e=-1`, and
:math:`{r_i}=0.01`. The thermal velocity is increased for ease of
visualization. Plot the simulated dispersion relation. You should see
another branch corresponding to the dispersion relation of waves
generated by the ions. Because they are more massive than electrons,
their frequency will be lower. These waves are called ion-acoustic waves
and the dispersion relation is given by:raw-latex:`\cite{Chen}`

.. math::

   \omega^2=\frac{\omega_{p,i}^2}{1+(\omega_{p,e}/k^2v_{\rm  th,e}^2)}.
               \label{IonDr}

 Plot the theoretical dispersion relations of
Eqs. \ `[LangDR] <#LangDR>`__ and `[IonDr] <#IonDr>`__ and compare them
with your simulations. Then decrease the thermal velocity of the
electrons, and determine how the dispersion relation changes.

*Problem 4: Buneman Instability*. What happens if there is a relative
velocity between two beams of different species? As in the
electron-electron two stream instability, an instability, known as the
Buneman instability, will be generated. It occurs if the drift velocity
between the electron beam and the ions exceeds the thermal velocity of
both species.:raw-latex:`\cite{Moreno2018ImpactOT}` Use the parameters
:math:`L=2\pi`, :math:`{N_t} = 5000`, :math:`\Delta t= 0.1`,
:math:`{N_g} =512`, :math:`N_e=10000`, :math:`N_i=10000`,
:math:`v0_{e}=1`, :math:`v0_{i}=0`, :math:`v_{\rm th,\,e}=0.004`,
:math:`v_{\rm th,\,i}=0`, :math:`r_e=-1`, :math:`r_i=0.001`. Plot the
phase space, and the spatial dependence of the electric potential and
electric field. Observe the time evolution and see how the plasma tries
to approach thermal equilibrium. You should first observe the formation
of electron phase-space holes with the characteristic bipolar waves of
the electric field propagating slower than in the electron-electron
two-stream instability.:raw-latex:`\cite{ionosfera}` You should also
observe how the electron velocity approximates a Maxwell-Boltzmann
distribution with a mean velocity approximately equal to the velocity of
the ion beam. This system requires a long time to reach thermal
equilibrium.

Next increase the charge-to-mass ratio of the ions. What happens to the
ion distribution compared to the first simulation? What would happen if
instead of ions, you consider a positron beam (:math:`r_i=1`)?

SG is grateful with Universidad EAFIT for its support. JH thanks the
support of University of Medellin, Colombia and Professor Jaime Araneda
of University of Concepción, Chile for his guidance in the first steps
of Particle simulations in Plasmas. JAV thanks the support of
ANID-Fondecyt under grant number 1190703.
