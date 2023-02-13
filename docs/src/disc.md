# All Functions with Hard-coded Discretizations

## Steady-State Heat Equation
- [`solve_T`](@ref)
    - Using 2nd order central differences, with linear extrapolation for ghost fluid cells.
    - Manually constructs sparse matrix, then solves it.
- [`solve_T_original`](@ref)
    - Using 2nd order central differences, with constant extrapolation for ghost fluid cells.
    - Manually constructs sparse matrix, then solves it.
- [`compute_frontvel_withT`](@ref)
    - Used together with [`vector_extrap_from_front`] to generate velocity field for advection
    - 1st order finite difference from interface cell into dried domain
- [`compute_frontvel_1`](@ref)
    - Used for generating arbitrary velocity fields--meant only for debugging
    - Used together with [`vector_extrap_from_front`] to generate velocity field for advection
    - 1st order finite difference from interface cell into dried domain
## Hyperbolic PDEs
- [`calc_Vd∇ϕ`](@ref)
    - Used in [`advect_ϕ!`](@ref) and [`advect_ϕ`](@ref)
    - Computes 1st order upwind finite differences
    - "Upwind": dictated by supplied velocity field
- [`calc_Nd∇v`](@ref) and [`calc_Nd∇v!`](@ref)
    - Used in [`vector_extrap_from_front`](@ref)
    - Unlike other cases, allocating and non-allocating are separate code
- [`𝒢`](@ref) , wrapped by [`𝒢_all`](@ref)
    - Used in [`reinitialize_ϕ!`](@ref)
    - Godunov's scheme, with first-order "upwind" discretization, for the norm of gradient
- [`calc_∇ϕ_1st`](@ref)
    - Not currently used
    - Compute 1st order upwind finite differences; roughly equivalent to Godunov's scheme (gave similar results, for uglier code)
    - "Upwind": towards interface
- [`calc_dϕdr_sdf`](@ref) and [`calc_dϕdz_sdf`](@ref)
    - Used inside [`update_ϕ_in_Γ!`](@ref), which is in turn inside [`reinitialize_ϕ!`](@ref)
    - Finite differences within interface region `Γ`; either first or second order, depending on availability
    - Implementation of Hartmann 2008 eq. 21, ignoring conditions which aid in coalescence
- [`identify_regions_RC`](@ref)
    - Second order finite differences, used inside [`update_ϕ_in_Γ!`](@ref) to estimate local curvature and identify cells which have more/less neighbors in interface
