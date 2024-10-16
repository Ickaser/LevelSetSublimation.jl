# All Functions with Hard-coded Discretizations

```@meta
CurrentModule = LevelSetSublimation
```

## Heat and Mass Transfer Equations
- [`solve_T`](@ref)
    - Using 2nd order central differences, with linear extrapolation for ghost fluid cells.
    - Manually constructs sparse matrix, then solves it.
    - Uses a 2nd-order discretization with linear extrapolation as defined in [gibouFourthOrderAccurate2005](@cite), section 2
- [`solve_p`](@ref)
## Velocity Calculations
- [`compute_frontvel_heat`](@ref)
    - Computes heat flux on Stefan interface
    - Uses 2nd-order finite difference with quadratic ghost cell extrapolation across Stefan interface to compute heat flux
    - Adapted from techniques explained in [gibouFourthOrderAccurate2005](@cite)
    - Generates a velocity field for level set advection, defined only on dry-layer side of interface; follow with `extrap_v_fastmarch` to extrapolate throughout domain
- [`compute_frontvel_mass`](@ref)
    - Computes heat flux on Stefan interface
    - Uses 2nd-order finite difference with quadratic ghost cell extrapolation across Stefan interface to compute heat flux
    - Adapted from techniques explained in [gibouFourthOrderAccurate2005](@cite)
    - Generates a velocity field for level set advection, defined only on dry-layer side of interface; follow with `extrap_v_fastmarch` to extrapolate throughout domain
- [`compute_frontvel_fixedspeed`](@ref)
    - Used for generating arbitrary velocity fields--meant only for debugging
    - Intended for use with velocity extrapolation; not currently implemented anywhere. 
    - 1st order finite difference from interface cell into dried domain
## Velocity Extrapolation
- [`extrap_v_fastmarch!`](@ref)
    - Uses 2nd-order upwind differences (towards the interface) in a fast-marching technique (which makes a single pass over all grid cells in sorted order) to extrapolate velocity away from the interface.
    - Follows [sethianFastMarchingMethods1999](@cite), section 5.1
## Reinitialization
- [`reinitialize_ϕ_HCR!`](@ref)
    - Reinitializes using a by-hand Explicit Euler in time and WENO in space, as described in [hartmannConstrainedReinitializationEquation2010](@cite)
## WENO Derivatives
- [`wenodiffs_local`](@ref)
    - Evaluates a single WENO derivative
    - Implemented as described in [hartmannAccuracyEfficiencyConstrained2009](@cite), 
- [`dϕdx_all_WENO`](@ref)
    - Used in [`𝒢_weno`](@ref) and therefore [`reinitialize_ϕ_HCR!`](@ref)
    - Also used in all velocity computation functions
- [`𝒢_weno`](@ref) , wrapped by `𝒢_weno_all`
    - Used in [`reinitialize_ϕ_HCR!`](@ref)
    - Godunov's scheme for the norm of gradient, using WENO derivatives
