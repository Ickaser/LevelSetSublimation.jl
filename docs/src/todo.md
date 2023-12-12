# Backlog

- Pull all level set geometry calculations into functions, write unit tests for all of those

- Build in some doctests
- Begin building test suite for every function, especially the edge cases
    - Velocity extrapolation
    - Geometry functions

- Reinit: steady state check for iteration
- Reinit happens at almost every time step. Speed up the check for determining if necessary?
- Docs for new reinit behavior in simulation
- Unify error region across all methods
- Use Della Rocca and Blanquart, 2014 for reinitialization at boundaries

- Begin reinvestigating performance, see if there are any easy gains to be had
- Examine current caching behavior
- Build a data dependency graph among functions, identify a better way to pass data around

- Strip Tf out of time integration, since it is not used

- Compile linsolve_testing_radial.jl, linsolve_testing_vertical.jl and pseudosteady_T_testing.jl into coherent unit tests; figure out a good place to keep the nice plots I generated
    - 1D temperatures
    - 1D pressures
    - Pseudosteady Tf
- Other cases to consider testing
    - 2D temperature, no ice

- Purge the many unnecessary lines of commented code

# Next TODO

- At boundaries, evaluate derivatives at interface instead of boundary (both T and p)

# Done (I think?)

- Write the "virtual thermocouple" function
- Formulate axisymmetric volume, surface area in terms of delta function integrals
- Build discrete delta integral
- Tests for discrete delta integral
- Incorporate discrete delta into Tf

- Make plots of LyoPRONTO comparison