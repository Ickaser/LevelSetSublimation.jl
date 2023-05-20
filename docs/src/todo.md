# To Do

- Use `u_from_ϕ_T` and similar functions everywhere, so that no state indices are hard-coded
- When building a "thermocouple" function for Tf, don't interpolate full system (just Tf at center), so that it's cheaper?
- Pull all level set geometry calculations into functions, write unit tests for all of those


- Begin adding variation in Tf: first, try with a full set of Tf variables, and see if it takes tiny time steps
    - Geometric: detecting where there's ice, what its height is, local surface area (mini contour)

- Build in some doctests
- Begin building test suite for every function, especially the edge cases
    - 1D temperatures
    - 2D temperature, no ice
    - 1D pressures, using flat/radial ice
    - Velocity extrapolation
    - Geometry functions

- Reinit: steady state check for iteration
- Reinit happens at almost every time step. Speed up the check for determining if necessary?
- Docs for new reinit behavior
- Unify error region across all methods