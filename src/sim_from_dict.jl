export sim_from_dict, sim_heatonly, sim_and_postprocess

# Callback functions ---------------
"""
    reinit_wrap(integ; verbose=false)

Thin wrapper to reinitialize the state of the level set function.

Uses a tolerance of 0.02 for error in norm of gradient (which should be 1)
 in the region B (band around the interface).
If `verbose=true`, logs the error in signed distance function before and after reinit, as well as current dried fraction
"""
function reinit_wrap(integ; verbose=false)
    dom = integ.p[1]
    ϕ = @views reshape(integ.u[iϕ(dom)], size(dom))
    verbose && (pre_err = sdf_err_L∞(ϕ, dom, region=:B))
    reinitialize_ϕ_HCR!(ϕ, dom, maxsteps=50, tol=0.02, err_reg=:B) 
    if verbose
        post_err = sdf_err_L∞(ϕ, dom, region=:B)
        dryfrac = 1 - compute_icevol_H(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        # @info "Reinit at t=$(integ.t)"
        @info "Reinit" integ.t pre_err post_err dryfrac spy(ϕ .< 0)
    end
end

function needs_reinit(u, t, integ)
    dom = integ.p[1]
    ϕ = reshape(u[iϕ(dom)], size(dom))
    # err = sdf_err_L1(ϕ, dom)
    B = identify_B(ϕ, dom)
    err = norm(𝒢_weno_all(ϕ, dom)[B].-1, Inf)
    tol = 0.02 # Based on Luddens 2015
    # @info "error from sdf" err-tol
    return err > tol
end

function store_Tf!(integ)
    Tf = @view integ.u[iTf(integ.p[1])]
    if integ.p[3]
        @info "callback" integ.t
        @time Tf_sol = pseudosteady_Tf(integ.u, integ.p[1], integ.p[2](integ.t), integ.p[4])
    else
        Tf_sol = pseudosteady_Tf(integ.u, integ.p[1], integ.p[2](integ.t), integ.p[4])
    end
    Tf .= Tf_sol
end

# Simulation functions -----------

"""
    sim_from_dict(fullconfig; tf=1e5, verbose=false)

Given a simulation configuration `fullconfig`, run a simulation.

Maximum simulation time (not CPU time, but simulation time) is specified by `tf`, which is in seconds; 1e5 is 270 hours. (No matter how much ice is left at that point, simulation will end.)
`verbose=true` will put out some info messages about simulation progress, i.e. at each reinitialization.

All passed parameters should have Unitful unit annotations, so they can be in whatever units are convenient;
inside this function, they will be converted to SI marks then have units stripped off (for easier numerical implementation).
The results will be in SI units, so times in seconds, temperature in Kelvin, pressure in Pa, etc.
Employ Unitful to do unit conversions after simulation if necessary.

`fullconfig` must have the following fields:
- `fillvol`, fill volume with Unitful units (e.g. `3u"mL"`)
- `vialsize`, a string (e.g. `"2R"`) giving vial size
- `paramsd`, a `Tuple{PhysicalProperties, TimeConstantProperties, TimeVaryingProperties}` 
    (i.e. a tuple of the three objects). See [`PhysicalProperties`](@ref), [`TimeConstantProperties`](@ref), [`TimeVaryingProperties`](@ref).
The following fields have default values and are therefore optional:
- `simgridsize`, a tuple giving number of grid points to use for simulation. Defaults to `(51, 51)`.
- `Tf0`, an initial ice temperature. Defaults to `Tsh(0)` if not provided  
- `Tvw0`, an initial glass temperature. Defaults to `Tf0`.
- `dudt_func`: defaults to [`dudt_heatmass!`](@ref); other options are [`dudt_heatmass_dae!`](@ref) and [`dudt_heatmass_implicit!`](@ref).
    Among these three, different problem formulations are used (explicit ODE with internal Newton solve, DAE, and implicit ODE).
- `init_prof`, a `Symbol` indicating a starting profile (from [`make_ϕ0`](@ref)

The three problem formulations have different advantages; in my testing, the DAE formulation and implicit formulation
tend to run faster and more closely reflect the problem structure,
but they run into instabilities when the sublimation front peels away from the wall,
whereas the explicit ODE formulation can jump over that point.
"""
function sim_from_dict(fullconfig; tf=1e5, verbose=false)

    # ------------------- Get simulation parameters

    @unpack paramsd = fullconfig
    # Default values for non-essential parameters
    init_prof = get(fullconfig, :init_prof, :flat) # Default to same ice & glass temperature if glass initial not given
    Tf0 = get(fullconfig, :Tf0, paramsd[3].Tsh(0u"s")) # Default to same ice & glass temperature if glass initial not given
    Tvw0 = get(fullconfig, :Tvw0, Tf0) # Default to same ice & glass temperature if glass initial not given

    # --------- Set up simulation domain, including grid size (defaults to 51x51)

    dom = Domain(fullconfig)
    # If no vial thickness defined, get it from vial size
    u0 = make_u0_ndim(init_prof, Tf0, Tvw0, dom)

    ϕ0 = @views reshape(u0[iϕ(dom)], size(dom))
    if verbose
        @info "Initializing ϕ"
    end
    # Make sure that the starting profile is very well-initialized
    # The chosen tolerance is designed to the error almost always seen in norm of the gradient
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=1000, tol=0.01, err_reg=:all) 

    sim_from_u0(u0, 0.0, fullconfig; tf=tf, verbose=verbose)
end

"""
    sim_from_u0(u0, t0, fullconfig; tf=1e5, verbose=false)

Wrapped by [`sim_from_dict`](@ref LevelSetSublimation.sim_from_dict); useful on its own if you want to start from partway through a simulation.
"""
function sim_from_u0(u0, t0, fullconfig; tf=1e5, verbose=false)

    @unpack paramsd = fullconfig
    dom = Domain(fullconfig)

    # ----- Nondimensionalize everything
    params_vary = params_nondim_setup(paramsd) # Covers the various physical parameters 

    # Cached array for using last pressure and Tf states as guess
    Tf0 = u0[iTf(dom)][1]
    Tf_last = fill(Tf0, dom.nr)


    dudt_func = get(fullconfig, :dudt_func, dudt_heatmass!) # Default to heat and mass transfer-based evolution
    # ---- Set up ODEProblem
    prob_pars = (dom, params_vary, verbose, Tf_last)
    tspan = (t0, tf)
    prob = ODEProblem(dudt_func, u0, tspan, prob_pars)

    # --- Set up reinitialization callback

    # After each time step, check if reinit is needed and carry out if necessary
    cb_reinit = DiscreteCallback(needs_reinit, x->reinit_wrap(x, verbose=verbose))

    # --- Set up simulation end callback

    # When the minimum value of ϕ is 0, front has disappeared
    cond_end(u, t, integ) = minimum(u) + min(dom.dz, dom.dr)/4 
    # ContinuousCallback gets thrown when `cond` evaluates to 0
    # `terminate!` ends the solve there
    cb_end = ContinuousCallback(cond_end, terminate!)

    # ------- Put together callbacks 
    # cbs = CallbackSet(cb_reinit, cb_end, cb_meas)
    cbs = CallbackSet(cb_end, cb_reinit)

    if verbose
        @info "Beginning solve"
    end

    # --- Solve
    if dudt_func == dudt_heatmass!
        # sol = solve(prob, SSPRK33(), dt=60, callback=cbs; ) # Fixed timestepping: 1 minute
        # sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping: default

        cb_store = DiscreteCallback((a,b,c)->true, store_Tf!, save_positions=(false,true))
        # @info "fixed"
        # sol = solve(prob, SSPRK33(), dt=60, callback=CallbackSet(cb_reinit, cb_end, cb_store); ) # Fixed timestepping: 1 minute
        sol = solve(prob, SSPRK43(), callback=CallbackSet(cb_reinit, cb_end, cb_store); ) # Adaptive timestepping: default
    elseif dudt_func == dudt_heatmass_dae!
        # Use a constant-mass-matrix representation with DAE, where Tf is algebraic
        massmat = Diagonal(vcat(ones(dom.nr*dom.nz), zeros(dom.nr), [1]))
        func = ODEFunction(dudt_heatmass_dae!, mass_matrix=massmat)
        prob = ODEProblem(func, u0, tspan, prob_pars)
        sol = solve(prob, FBDF(); callback=cbs)
    elseif dudt_func == dudt_heatmass_implicit!
        sol = solve(prob, Rodas4P(), callback=cbs; ) # Adaptive timestepping: default
    else
        sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping: default
    end
    return @strdict sol dom
end

function sim_and_postprocess(config; kwargs...)
    @time res = sim_from_dict(config; kwargs...)
    rpos = [0,   0  , 0, 0.5, 1]
    zpos = [0.5, 0.25, 0, 0  , 0]
    @time Tf_sol = virtual_thermocouple(rpos, zpos, res, config)
    return @strdict res Tf_sol config
end


