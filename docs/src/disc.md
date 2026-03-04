# All Functions with Hard-coded Discretizations

```@meta
CurrentModule = LevelSetSublimation
```

## Heat and Mass Transfer Equations
- [`solve_T`](@ref)
    - Using 2nd order central differences, with linear extrapolation for ghost fluid cells across the level set interface.
    - Manually constructs sparse matrix, then solves it.
    - Uses a 2nd-order discretization with linear extrapolation as defined in [gibouFourthOrderAccurate2005](@cite), section 2
    - All the possible stencils for this discretization were generated in a Jupyter notebook, visible in this code's GitHub repo at `/docs/gfm_extrap.ipynb`. If links don't break, that can be viewed online [here](https://nbviewer.org/github/ickaser/LevelSetSublimation.jl/blob/main/docs/src/gfm_extrap.ipynb).
- [`solve_p`](@ref)
    - Since the mass conductivity in the equation for $p$ may itself depend on $p$, this function guesses values of $b$, solves for $p$, then iterates a few more times.
- [`solve_p_given_b`](@ref)
    - This does the actual linear solve for `solve_p`. Like `solve_T`, it uses a 2nd-order discretization as in [gibouFourthOrderAccurate2005](@cite), section 2. This becomes slightly more complicated though because the mass conductivity $b$ is not spatially constant.
    - All the possible stencils for this discretization were generated in a Jupyter notebook, visible in this code's GitHub repo at `/docs/gfm_extrap_vark.ipynb`. If links don't break, that can be viewed online [here](https://nbviewer.org/github/ickaser/LevelSetSublimation.jl/blob/main/docs/src/gfm_extrap_vark.ipynb).
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
