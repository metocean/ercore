Introduction
============

ERCore is a Lagrangian particle tracking model for simulating the
dispersion and fate of material in the oceanic or atmospheric
environment. It solves equations describing the movement and
modification of discrete particles in response to the environment. The
ambient flow of water or wind advects the particles with spatially and
temporally varying fields. Turbulent diffusion is included with a random
walk approach and the presence of a boundary such as a shoreline or the
air-sea interface can be specified. A simulated particle can represent
any material with additional physical, chemical or biological processes
added as required to simulate active modifications that may occur as it
moves within the water or air. Examples of active particles include
sediment, oil, plankton and even persons lost at sea.

This document describes the mathematical and numerical formulation of
the model. Also included are descriptions of the various extensions that
can be applied to the particles to allow them to simulate different
materials.

Numerical formulation
=====================

Equations of motion
-------------------

The coordinate system is defined with $z$ positive upwards. The nominal
0 level for $z$ is Mean Sea Level, however any consistent vertical
reference can be used. The change of position $x,y,z$ of a particle at
time $t$ is governed by:$$\begin{aligned}
\frac{dx}{dt} = U(x,y,z,t)+\tilde{u}+u_p \\
\frac{dy}{dt} = V(x,y,z,t)+\tilde{v}+v_p \\ 
\frac{dz}{dt} = W(x,y,z,t)+\tilde{w}+w_p \end{aligned}$$ where $U,V,W$
is the ambient flow and $\tilde{u},\tilde{v},\tilde{w}$ are a turbulent
flow component. $u_p,v_p,w_p$ are any additional motions that are
properties of the particle’s interaction with the surrounding water or
air such as sinking, floating, drifting or swimming.

Each particle has a weighting and age attribute attached to it. The
weight value can represent any measure of quantity depending on the
application, for example mass, volume or number. The generic equation
describing the rate of change of weighting, $N$ is:
$$\frac{dN}{dt}=f(P,U,V,W,S,T,D,H)$$ where $P$ are properties of the
particle itself and $S,T,D,H$ are salinity, temperature, density and
humidity of water or air. Note that for a given process and medium only
some of the ambient properties will be relevant. The age is a simple
measure of time since a particle was released. Any number of other
attributes can also be attached to a particle which describe its
physical, chemical or biological state.

![Schematic of physical processes acting on a particle in the
model.](model.png){width="\textwidth"}

Turbulent diffusion
-------------------

The turbulent component of the velocity is given by: $$\begin{aligned}
\int_{t}^{t+\Delta t} \tilde{u} dt = \sqrt{6K_h \Delta t} H(-1,1) \\
\int_{t}^{t+\Delta t} \tilde{v} dt = \sqrt{6K_h \Delta t} H(-1,1) \\
\int_{t}^{t+\Delta t} \tilde{w} dt = \sqrt{6K_v \Delta t} H(-1,1)\end{aligned}$$
where $K_h$ and $K_v$ are horizontal and vertical eddy diffusion
coefficients and $H$ is a number sampled from a uniform distribution
between $-1$ and $1$. The eddy diffusion coefficients can be provided as
known values, either constant or as a spatially (and temporally) varying
field.

Coordinate system
-----------------

The model solves all equations in standard SI units. However the
horizontal coordinates can also be configured as geographical longitude
and latitude. In this case a map factor is applied to any horizontal
differences as: $$\begin{aligned}
d\lambda = \frac{360dx}{C*cos{\theta}} \\
d\theta = \frac{360dy}{C}\end{aligned}$$ where $\lambda$ and $\theta$
are the longitude and latitude respectively. $C$ is the circumference of
the earth.

Numerical solution
------------------

The equations of motion for each particle are integrated with a constant
time step, $\Delta t$. A Runge-Kutta scheme is used for the ambient flow
component, with default order 4 (which can be reduced for computational
speed by user configuration). The turbulent and particle components of
motion are added with a simple first order integration.

Forcing fields
--------------

The ambient flow $U,V,W$ is provided as input data which can be constant
or as spatially and temporally varying fields. For variable fields, the
value at the particle position for a given time is determined by simple
multilinear interpolation. Any number of additional fields which can
influence the particle can also be prescribed and are likewise
interpolated to provide the local ambient conditions for the particle at
each timestep.

Boundaries
----------

Boundaries can be specified in 4 ways:

1.  A shoreline vector

2.  A spatially varying field (such as ocean depth)

3.  A constant level

4.  A bounding box

In each case an intersection algorithm is applied at each time-step to
determine whether the particle crosses the boundary. Firstly the
particle motion vector for the time-step is calculated, then an
intersect test is made between that vector and the boundaries.

Shorelines are defined as line segments, each of which is tested for
intersection with the motion vector in the horizontal plane. Shorelines
are assumed to exist at a z-level of zero, and any horizontal
intersection will relocate the particle at that level.

For variable boundaries such as the sea bottom, the z-coordinate at the
start and end of the motion vector are first interpolated. The straight
line between these two z-coordinates is then used for the intersection
test with the motion vector in the vertical plane; the horizontal
position of an intersect is then the fractional displacement in the
horizontal direction.

If an intersection occurs, the particle is placed at the intersection
point. The particle’s subsequent action at the following time-step
depends on whether the boundary is sticky. If the boundary is non-sticky
the particle is free to move away from the boundary if its vector of
motion carries it away, otherwise it remains at its intersection
location. In the case of a bounding box, the particle is removed from
the simulation.

![Illustration of boundary crossing tests for shoreline (top) and
variable surface (bottom).](intersect.png){width="\textwidth"}

Particle release and tracking
-----------------------------

Particle releases are specified to be a particular location and the
following options are available:

1.  Instantaneaous release of total material $M_{tot}$: All $n$
    particles added at a given time. $$N=\frac{N_tot}{n}$$

2.  Continuous release, with flux of material $F$: Particles are added
    continuously, at a rate of $n_{rel}$ every timestep
    $\Delta t_{rel}$. $$N=\frac{F}{\Delta t_{rel} n_{rel}}$$

3.  Staged release of total material $N_{tot}$ added over time.
    $$N=\frac{N_{tot}}{\Delta t_{rel}}$$

Multiple releases can be specified and subsequently tracked. A status
code on each particle flags its state, for example whether it is
releases and mobile, stuck to a boundary or has been removed. The
particle location and statuses are output at user specified intervals.

Particle Processes
==================

Biota
-----

Processes currently included that affect biota are:

1.  Mortality in response to age, temperature or salinity

2.  Spawning to become a subsequent life stage

3.  Diurnal vertical migration

Often, the weight property will represent number of organisms that each
particle represents.

Mortality is modelled with a simple die-off rate: $$\frac{dW}{dt}=-rW$$
where $r$ is the mortality rate. The mortality rate, can either be
specified as a constant value or as a function of temperature and/or
salinity: $$\begin{aligned}
r=\frac{T-T_L^{high}+T_{tol}}{T_{tol}} \in[0,1] \\
r=\frac{T_L^{low}+T_{tol}-T}{T_{tol}} \in[0,1]\end{aligned}$$ where
$T_L$ is the lethal temperature and $T_{tol}$ is a tolerance which
describes the temperature range over which effects start to be felt.
Exactly analogous equations exist for salinity. When the particle weight
drops below a prescribed number, the particle is removed from the
simulation.

Diurnal vertical migration is determined by: $$\begin{aligned}
w_p=w_m(t) \\
z \in [-z_{day},-z_{night}]\end{aligned}$$ where $w_m$ is a vertical
migration constant which is positive (up) after the sun sets, and
negative down after sun rise.

Spawning occurs when biota either mature to their next life stage or
release offspring. In the first case the relevant class of the biota
simply changes. In the latter a new particle is created at the existing
position, which its weight determined by a spawning ratio that describes
the

Sediment
--------

Simulation of sediment in water or particles in air can incorporate a
density based calculation of fall velocity rather than a constant
specified value. $$w_p=$$

Oil Spill
---------

Configuration
=============

Tests
=====

Comparison with GNOME
---------------------

GNOME is the NOAA oil spill model
(<http://response.restoration.noaa.gov/gnome>), which follows the same
essential methodology as ERCore. GNOME simulates particle motion in a
temporally and spatially variable horizontal field. To verify that
ERCore and GNOME produce the same trajectories in the absence of
turbulent diffusion, both models were run with a single particle in an
identical two-dimensional flow field. The figure below shows that both
models produce identical particle tracks in a complex flow field.

![Single particle trajectories from GNOME and ERCore for an identical
flow field in the absence of
diffusion](ercoregnome.png){width="\textwidth"}
