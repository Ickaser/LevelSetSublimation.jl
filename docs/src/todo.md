# To Do

- Create a function `u_from_ϕ_T` so that no state indices are hard-coded
- When building a "thermocouple" function for Tf, don't interpolate full system (just Tf at center), so that it's cheaper?


- Begin adding variation in Tf: first, try with a full set of Tf variables, and see if it takes tiny time steps
    - Geometric: detecting where there's ice, what its height is, local surface area (mini contour)

- Build in some doctests
- Begin building test suite for every function, especially the edge cases
    - Start from easy functions and early-defined functions
    - Consider validation cases from early notebooks

- Extrapolation for WENO: try removing the boundary direction specialization where that was used in choice of east/west derivatives (left boundary can use all extrapolated cells, if velocity dictates)

- G for WENO: improve performance
- Reinit conditions: play with tolerances
- Reinit: steady state check for iteration
- Velocity extrapolation stil needs behavior checking

- Reinit happens at almost every time step. Speed up the need to check it?
- Docs for new reinit
- Unify error region across all methods