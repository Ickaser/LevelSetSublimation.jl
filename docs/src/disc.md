# All Functions with Hard-coded Discretizations

## Steady-State Heat Equation
- [`solve_T`](@ref)
    - Using 2nd order central differences, with linear extrapolation for ghost fluid cells.
    - Manually constructs sparse matrix, then solves it.
    - Uses a 2nd-order discretization with linear extrapolation as defined in [gibouFourthOrderAccurate2005](@cite), section 2
- `solve_T_original`
    - Using 2nd order central differences, with constant extrapolation for ghost fluid cells.
    - Manually constructs sparse matrix, then solves it.
    - Equivalent to a 2nd-order discretization with constant extrapolation as defined in [gibouFourthOrderAccurate2005](@cite)
    - (Globally, 1st order accurate.)
- [`compute_frontvel_withT`](@ref)
    - Computes heat flux on Stefan interface
    - Uses 2nd-order finite difference with quadratic ghost cell extrapolation across Stefan interface to compute heat flux
    - Used by [`extrap_v_pde`](@ref) or [`extrap_v_fastmarch`](@ref) to generate velocity field for level set advection
    - Adapted from techniques explained in [gibouFourthOrderAccurate2005](@cite)
## Velocity Extrapolation
- [`extrap_v_fastmarch`](@ref)
    - Uses 2nd-order upwind differences (towards the interface) in a fast-marching technique (which makes a single pass over all grid cells in sorted order) to extrapolate velocity away from the interface.
    - Follows [sethianFastMarchingMethods1999](@cite), section 5.1
- [`compute_frontvel_1`](@ref)
    - Used for generating arbitrary velocity fields--meant only for debugging
    - Intended for use with velocity extrapolation; not currently implemented anywhere. 
    - 1st order finite difference from interface cell into dried domain
## Reinitialization
- [`𝒢_weno`](@ref) , wrapped by `𝒢_weno_all`
    - Used in [`reinitialize_ϕ_HCR!`](@ref)
    - Godunov's scheme for the norm of gradient, using WENO derivatives
- [`𝒢_1st`](@ref) , wrapped by `𝒢_1st_all`
    - Used in [`reinitialize_ϕ!`](@ref)
    - Godunov's scheme for the norm of gradient using first-order upwind derivatives
## WENO Derivatives
- [`wenodiffs_local`](@ref)
    - Implemented as described in [hartmannConstrainedReinitializationEquation2010](@cite)
    - Used in [`𝒢_weno`](@ref) and therefore [`reinitialize_ϕ_HCR!`](@ref)
    - Used in [`]

## Obsolete Methods
- [`calc_Vd∇ϕ`](@ref)
    - Used in [`advect_ϕ!`](@ref) and [`advect_ϕ`](@ref)
    - Computes 1st order upwind finite differences
    - "Upwind": dictated by supplied velocity field
- [`extrap_v_pde`](@ref)
    - Less accurate and more expensive than [`extrap_v_fastmarch`](@ref)
    - Uses 1st-order upwind differences, in theory towards the interface, and advances a PDE in pseudo-time to generate the desired extrapolation.
    - Drawn from [osherLevelSetMethods2003](@cite)
- [`calc_Nd∇v`](@ref) and [`calc_Nd∇v!`](@ref)
    - Used in [`extrap_v_pde`](@ref)
    - Unlike other cases, allocating and non-allocating are separate code
- [`calc_∇ϕ_1st`](@ref)
    - Compute 1st order upwind finite differences, upwind towards interface; replaced by 
- [`calc_dϕdr_sdf`](@ref) and [`calc_dϕdz_sdf`](@ref)
    - Used inside [`update_ϕ_in_Γ!`](@ref), which is in turn inside [`reinitialize_ϕ!`](@ref)
    - Finite differences within interface region `Γ`; either first or second order, depending on availability
    - Implementation of [hartmannDifferentialEquationBased2008](@cite) eq. 21, ignoring conditions which aid in coalescence
- [`identify_regions_RC`](@ref)
    - Second order finite differences, used inside [`update_ϕ_in_Γ!`](@ref) to estimate local curvature and identify cells which have more/less neighbors in interface
