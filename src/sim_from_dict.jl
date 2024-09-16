export sim_from_dict, sim_heatonly, sim_and_postprocess

"""
    reinit_wrap(integ; verbose=false)

Thin wrapper to reinitialize the state of the level set function.

Uses a tolerance of 0.02 for error in norm of gradient (which should be 1)
 in the region B (band around the interface).
If `verbose=true`, logs the error in signed distance function before and after reinit, as well as current dried fraction
"""
function reinit_wrap(integ; verbose=false)
    dom = integ.p[1]
    ϕ = ϕ_T_from_u_view(integ.u, dom)[1]
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
    ϕ = ϕ_T_from_u(u, dom)[1]
    # err = sdf_err_L1(ϕ, dom)
    B = identify_B(ϕ, dom)
    err = norm(𝒢_weno_all(ϕ, dom)[B].-1, Inf)
    tol = 0.02 # Based on Luddens 2015
    # @info "error from sdf" err-tol
    return err > tol
end

function input_measurements!(params, t::Number, controls)
    for key in keys(controls)
        params[key] = controls[key](t)
    end
end

"""
    sim_from_dict(fullconfig; tf=1e5, verbose=false)

*TODO out of date!!!!*
Given a simulation configuration `fullconfig`, run a simulation.

Maximum simulation time (not CPU time, but simulation time) is specified by `tf`, which is in seconds; 1e5 is 270 hours. (No matter how much ice is left at that point, simulation will end.)
`verbose=true` will put out some info messages about simulation progress, i.e. at each reinitialization.

All passed parameters should have Unitful unit annotations, so they can be in whatever units are convenient;
inside this function, they will be converted to SI marks then have units stripped off (for easier numerical implementation).
The results will be in SI units, so times in seconds, temperature in Kelvin, pressure in Pa, etc.
Employ Unitful to do unit conversions after simulation if necessary.

`fullconfig` should have the following fields:
- `init_prof`, types listed for [`make_ϕ0`](@ref)
- `fillvol`, fill volume with Unitful units (e.g. `3u"mL"`)
- `vialsize`, a string (e.g. `"2R"`) giving vial size
- `simgridsize`, a tuple/arraylike giving number of grid points to use for simulation. Defaults to `(51, 51)`.
- `Tf0`, an initial ice temperature with Unitful units 
- `Tvw0`, an initial glass temperature (if the same as Tf0, can leave this out)
- `cparams`, which in turn has fields with Unitful units
    - `Kvwf`, 
    - `Kshf` : heat transfer coefficients shelf
    # - `QRFd` : volumetric heating in cake 
    - `kd`: thermal conductivity of cake
    - `m_cp_gl` total thermal mass of glass, relevant to heating/cooling of glass wall
    - `kf`: thermal conductivity of ice
    - `ρf`: density of ice
    - `Cpf`: heat capacity of ice
    - `ΔH` : heat of sublimation of ice
    - `ϵ` : porosity of porous medium
    - `l` : dusty gas model constant: characteristic length for Knudsen diffusion
    - `κ` : dusty gas model constant: length^2 corresponding loosely to Darcy's Law permeability
    - `Mw`: molecular weight of species (water), with appropriate units
    - `μ` : dynamic viscosity of species (water), with appropriate units
- `dudt_func`: defaults to `dudt_heatmass!`, but for other cases (e.g. ignoring mass transfer) can replace this

If you are getting a warning about instability, it can sometimes be fixed by tinkering with the reinitialization behavior.


"""
function sim_from_dict(fullconfig; tf=1e5, verbose=false)

    # ------------------- Get simulation parameters

    @unpack init_prof, Tf0, = fullconfig
    # Default values for non-essential parameters
    Tvw0 = get(fullconfig, :Tvw0, Tf0) # Default to same ice & glass temperature if glass initial not given

    # --------- Set up simulation domain, including grid size (defaults to 51x51)

    dom = Domain(fullconfig)
    # If no vial thickness defined, get it from vial size
    u0 = make_u0_ndim(init_prof, Tf0, Tvw0, dom)

    ϕ0 = ϕ_T_from_u_view(u0, dom)[1]
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

Wrapped by [`sim_from_config`](@ref); useful on its own if you want to start from partway through a simulation.
"""
function sim_from_u0(u0, t0, fullconfig; tf=1e5, verbose=false)

    @unpack paramsd = fullconfig
    dom = Domain(fullconfig)

    # ----- Nondimensionalize everything
    params_vary = params_nondim_setup(paramsd) # Covers the various physical parameters 

    # Cached array for using last pressure and Tf states as guess
    # This gets ignored for heat-only simulation, but shouldn't cause any problems
    Tf0 = ϕ_T_from_u(u0, dom)[2][1]
    p_sub = calc_psub(Tf0)
    p_last = fill(p_sub, size(dom))
    Tf_last = fill(Tf0, dom.nr)


    dudt_func = get(fullconfig, :dudt_func, dudt_heatmass!) # Default to heat and mass transfer-based evolution
    # ---- Set up ODEProblem
    prob_pars = (dom, params_vary, p_last, Tf_last, verbose)
    tspan = (0, tf)
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
        function store_Tf!(integrator)
            ϕ, Tf, Tvw = ϕ_T_from_u_view(integrator.u, integrator.p[1])
            params = integrator.p[2](integrator.t)
            verbose && @info "callback" integrator.t
            @time Tf_sol = pseudosteady_Tf(integrator.u, integrator.p[1], params, integrator.p[4])
            Tf .= Tf_sol
        end

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


