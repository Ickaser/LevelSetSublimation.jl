# Backlog

- Pull all level set geometry calculations into functions, write unit tests for all of those

- Build in some doctests
- Begin building test suite for every function, especially the edge cases
    - 1D temperatures
    - 2D temperature, no ice
    - 1D pressures, using flat/radial ice
    - Velocity extrapolation
    - Geometry functions

- Reinit: steady state check for iteration
- Reinit happens at almost every time step. Speed up the check for determining if necessary?
- Docs for new reinit behavior in simulation
- Unify error region across all methods
- Use Della Rocca and Blanquart, 2014 for reinitialization at boundaries

# Next TODO



- Change process conditions (Tsh, Pch, etc.) into a time-dependent function, rather than a callback (currently holds back the time stepping)
- Put RampedVariables in as time-dependent part of ODE function

- Improve speed of pseudosteady calculation

# Things to clean up later
- Passing structure inside time derivatives and simulation "parameters"

# Done (I think?)

- Formulate axisymmetric volume, surface area in terms of delta function integrals
- Build discrete delta integral
- Tests for discrete delta integral
- Incorporate discrete delta into Tf

- Make plots of LyoPRONTO comparison