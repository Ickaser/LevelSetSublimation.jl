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
        dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        # @info "Reinit at t=$(integ.t)"
        @info "Reinit" integ.t pre_err post_err dryfrac spy(ϕ .< 0)
    end
end

# """
# next_reinit_time(integ)

# Compute the next time reinitialization should be necessary, given the current integrator state.
# Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
# """
# function next_reinit_time(integ; verbose=false)
#     dom = integ.p[1]
#     du = similar(integ.u)
#     dudt_heatmass!(du, integ.u, integ.p, integ.t)

#     # The main region of concern is the frozen region near interface
#     # Find the largest value of dϕdt in that region
#     ϕ, Tf, Tw = ϕ_T_from_u(integ.u, dom)
#     dϕ, dTfdt, dTwdt = ϕ_T_from_u(du, dom)
#     B = identify_B(ϕ, dom)
#     B⁻ = B .& (ϕ .<= 0) # BitArray: Frozen region, near interface
#     if sum(B⁻) > 0 # FOund a frozen region
#         max_dϕdt = maximum(abs.(dϕ[B⁻]))
#     else # No frozen region: use whole domain
#         max_dϕdt = maximum(abs.(dϕ))
#     end
#     @unpack p_ch = integ.p[2]
#     # if max_dϕdt == 0 # If no sublimation is happening...
#     if calc_psub(Tf) < p_ch
#         # Guess when it will start
#         # Find temperature such that p_sub = p_ch
#         Tstart = Tf
#         for i in 1:5 # 5 steps of Newton iteration
#             pg = calc_psub(Tstart)
#             dpdT = pg - calc_psub(Tstart - 1)
#             dT = -(pg-p_ch)/dpdT
#             Tstart += dT
#         end
#         dt = abs((Tstart - Tf)/dTfdt) * 1.05 # Overshoot slightly, to get into real drying behavior before next check time
#         # @info "Sublimation stopped, estimating reinit." calc_psub(Tstart) calc_psub(Tf) p_ch dt dTfdt
#         return integ.t + dt
#     end

#     # Reinit next when interface should have moved across half of band around interface 
#     domfrac = 0.1 # Should reinit roughly 10-15 times during simulation
#     minlen = min(domfrac*dom.rmax , domfrac*dom.zmax)  
#     dt = minlen / max_dϕdt
#     dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
#     if dryfrac < domfrac || integ.t < 1000
#         dt = clamp(dt, 1, 100)
#     end
#     # @info "Reinit at t=$(integ.t), dt=$dt" minlen extrema(dϕ[B⁻])#, next at t=$(integ.t+dt)" 
#     if verbose
#         reinit_err = sdf_err_L∞(ϕ, dom)
#         # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
#         @info "Reinit at t=$(integ.t), dt=$dt" dryfrac Tf Tw reinit_err #, next at t=$(integ.t+dt)" 
#     end

#     return integ.t + dt
# end

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
- `Tw0`, an initial glass temperature (if the same as Tf0, can leave this out)
- `controls`, which has following fields (either scalar or array, with same length as `t_samp`:
    - `t_samp`, sampled measurement times. Needed only if other measurements are given during time
    - `Tsh`, shelf temperature: either a scalar (constant for full time span) or an array at specified time, in which case implemented via callback
    - `Q_gl_RF`, glass RF heating. Scalar or array, like Tsh
    - `Q_ic`, ice RF heating. 
    - `p_ch` : pressure at top of cake
- `cparams`, which in turn has fields with Unitful units
    - `Kw`, 
    - `Kv` : heat transfer coefficients shelf
    - `Q_ck` : volumetric heating in cake 
    - `k`: thermal conductivity of cake
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

During simulation, at each value of `t_samp`, the values of any `controls` which are arrays will be added to an internal dict called `params`.

If you pass in an array of values for multiple of `Tsh`, `Q_gl_RF`, or others, they must all have the same length as `t_samp`.

If you are getting a warning about instability, it can sometimes be fixed by tinkering with the reinitialization behavior.


"""
function sim_from_dict(fullconfig; tf=1e5, verbose=false)

    # ------------------- Get simulation parameters

    @unpack cparams, init_prof, Tf0, controls, vialsize, fillvol = fullconfig

    # Default values for non-essential parameters
    Tw0 = get(fullconfig, :Tw0, Tf0) # Default to same ice & glass temperature if glass initial not given
    dudt_func = get(fullconfig, :dudt_func, dudt_heatmass!) # Default to heat and mass transfer-based evolution

    # --------- Set up simulation domain, including grid size (defaults to 51x51)

    dom = Domain(fullconfig)
    # If no vial thickness defined, get it from vial size
    if !haskey(cparams, :vial_thick) 
        cparams[:vial_thick] =  get_vial_thickness(vialsize)
    end

    # ----- Nondimensionalize everything

    params, ncontrols = params_nondim_setup(cparams, controls) # Covers the various physical parameters 

    u0 = make_u0_ndim(init_prof, Tf0, Tw0, dom)

    ϕ0 = ϕ_T_from_u_view(u0, dom)[1]
    if verbose
        @info "Initializing ϕ"
    end
    # Make sure that the starting profile is very well-initialized
    # The chosen tolerance is designed to the error almost always seen in norm of the gradient
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=1000, tol=0.01, err_reg=:all) 

    # Cached array for using last pressure and Tf states as guess
    # This gets ignored for heat-only simulation, but shouldn't cause any problems
    p_sub = calc_psub(ustrip(u"K", Tf0))
    p_last = fill(p_sub, size(dom))
    Tf_last = fill(ustrip(u"K", Tf0), dom.nr)


    # ---- Set up ODEProblem
    prob_pars = (dom, params, p_last, Tf_last, ncontrols, verbose)
    tspan = (0, tf)
    prob = ODEProblem(dudt_func, u0, tspan, prob_pars)

    # --- Set up reinitialization callback

    # cb_reinit = IterativeCallback(x->next_reinit_time(x, verbose=verbose), reinit_wrap,  initial_affect = true)
    # After each time step, check if reinit is needed and carry out if necessary
    cb_reinit = DiscreteCallback(needs_reinit, x->reinit_wrap(x, verbose=verbose))

    # --- Set up simulation end callback

    # When the minimum value of ϕ is 0, front has disappeared
    cond_end(u, t, integ) = minimum(u) + min(dom.dz, dom.dr)/4 
    # ContinuousCallback gets thrown when `cond` evaluates to 0
    # `terminate!` ends the solve there
    cb_end = ContinuousCallback(cond_end, terminate!)

    # # Check after each time step, rather than seeking exact end time
    # cond_end(u, t, integ) = all(u .>= 0)
    # cb_end = DiscreteCallback(cond_end, terminate!)

    # ----- Set up  measurement callback
    # meas_affect!(integ) = input_measurements!(integ, meas_keys, ncontrols)
    # cb_meas = PresetTimeCallback(ncontrols[:t_samp], meas_affect!, filter_tstops=true)

    # ------- Put together callbacks 
    # cbs = CallbackSet(cb_reinit, cb_end, cb_meas)
    cbs = CallbackSet(cb_reinit, cb_end)

    if verbose
        @info "Beginning solve"
    end
    # --- Solve
    if dudt_func == dudt_heatmass!
        # sol = solve(prob, SSPRK33(), dt=60, callback=cbs; ) # Fixed timestepping: 1 minute
        # sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping: default
        function store_Tf!(integrator)
            ϕ, Tf, Tw = ϕ_T_from_u_view(integrator.u, integrator.p[1])
            verbose && @info "callback" integrator.t
            @time Tf_sol = pseudosteady_Tf(integrator.u, integrator.p[1], integrator.p[2], integrator.p[4])
            Tf .= Tf_sol
        end

        cb_store = DiscreteCallback((a,b,c)->true, store_Tf!, save_positions=(false,true))
        # @info "fixed"
        # sol = solve(prob, SSPRK33(), dt=60, callback=CallbackSet(cb_reinit, cb_end, cb_store); ) # Fixed timestepping: 1 minute
        sol = solve(prob, SSPRK43(), callback=CallbackSet(cb_reinit, cb_end, cb_store); ) # Adaptive timestepping: default
    elseif dudt_func == dudt_heatmass_dae! || dudt_func == dudt_heatmass_dae_or_newton!
        # Use a constant-mass-matrix representation with DAE, where Tf is algebraic
        massmat = Diagonal(vcat(ones(length(ϕ0)), zeros(dom.nr), [1]))
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

"""
Flat copy of sim_from_dict as of 2024-02-07, except it takes u0 as argument, allowing starting from midway through cycle.
"""
function sim_from_u0(u0, t0, fullconfig; tf=1e5, verbose=false)
    @unpack cparams, init_prof, controls, vialsize, fillvol = fullconfig
    dudt_func = get(fullconfig, :dudt_func, dudt_heatmass!) # Default to heat and mass transfer-based evolution
    dom = Domain(fullconfig)
    if !haskey(cparams, :vial_thick) 
        cparams[:vial_thick] =  get_vial_thickness(vialsize)
    end
    params, ncontrols = params_nondim_setup(cparams, controls) # Covers the various physical parameters 
    ϕ0, Tf0, Tvw0 = ϕ_T_from_u_view(u0, dom)
    if verbose
        @info "Initializing ϕ"
    end
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=1000, tol=0.01, err_reg=:all) 
    # p_sub = calc_psub(ustrip(u"K", Tf0))
    p_sub = calc_psub(Tf0[1])
    p_last = fill(p_sub, size(dom))
    # Tf_last = fill(ustrip(u"K", Tf0), dom.nr)
    Tf_last = copy(Tf0)
    prob_pars = (dom, params, p_last, Tf_last, ncontrols, verbose)
    tspan = (t0, tf)
    prob = ODEProblem(dudt_func, u0, tspan, prob_pars)
    cb_reinit = DiscreteCallback(needs_reinit, x->reinit_wrap(x, verbose=verbose))
    cond_end(u, t, integ) = minimum(u) 
    cb_end = ContinuousCallback(cond_end, terminate!)
    cbs = CallbackSet(cb_reinit, cb_end)
    if verbose
        @info "Beginning solve" u0 t0
    end
    if dudt_func == dudt_heatmass!
        function store_Tf!(integrator)
            ϕ, Tf, Tw = ϕ_T_from_u_view(integrator.u)
            @info "callback"
            @time Tf_sol = pseudosteady_Tf(integrator.u, integrator.p[1], integrator.p[2], integrator.p[4])
            Tf .= Tf_sol
        end

        cb_store = DiscreteCallback((a,b,c)->true, store_Tf!, store_positions=(false,true))
        # sol = solve(prob, SSPRK33(), dt=60, callback=cbs; ) # Fixed timestepping: 1 minute
        sol = solve(prob, SSPRK43(), callback=CallbackSet(cb_reinit, cb_end, cb_store); ) # Adaptive timestepping: default
    elseif dudt_func == dudt_heatmass_dae! || dudt_func == dudt_heatmass_dae_or_newton!
        # Use a constant-mass-matrix representation with DAE, where Tf is algebraic
        massmat = Diagonal(vcat(ones(length(ϕ0)), zeros(dom.nr), [1]))
        func = ODEFunction(dudt_func, mass_matrix=massmat)
        prob = ODEProblem(func, u0, tspan, prob_pars)
        sol = solve(prob, Rodas4P2(); callback=cbs)
        # sol = solve(prob, FBDF(); callback=cbs)
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



# """
# next_reinit_time_heatonly(integ)

# Compute the next time reinitialization should be necessary, given the current integrator state.
# Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
# """
# function next_reinit_time_heatonly(integ; verbose=false)
#     dom = integ.p[1]
#     du = similar(integ.u)
#     dudt_heatonly!(du, integ.u, integ.p, integ.t)

#     # The main region of concern is the frozen region near interface
#     # Find the largest value of dϕdt in that region
#     ϕ, Tf, Tw = ϕ_T_from_u(integ.u, dom)
#     dϕ, dTfdt, dTwdt = ϕ_T_from_u(du, dom)
#     B = identify_B(ϕ, dom)
#     B⁻ = B .& (ϕ .<= 0)
#     if sum(B⁻) > 0
#         max_dϕdt = maximum(abs.(dϕ[B⁻]))
#     else
#         max_dϕdt = maximum(abs.(dϕ))
#     end
#     @unpack p_ch = integ.p[2]
#     if max_dϕdt == 0 # If no sublimation is happening...
#         @info "Sublimation stopped, no way to estimate next reinit time so guessing dt=1.0 ." 
#         return integ.t + 1.0
#     end

#     # Reinit next when interface should have moved across half of band around interface 
#     domfrac = 0.1 # Should reinit roughly 10-15 times during simulation
#     minlen = min(domfrac*dom.rmax , domfrac*dom.zmax)  # Also: at early times, do more often
#     dt = minlen / max_dϕdt
#     dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
#     if dryfrac < domfrac #|| integ.t < 1000
#         # dt = clamp(dt, 1, )
#         # integ.t / dryfrac * 0.1 / 2
#         integ.t / dryfrac * 0.01 / 2
#         dt = max(1.0, integ.t / (dryfrac-1e-4) * 0.05) # Expect another .05 dryfrac
#     end
#     # @info "Reinit at t=$(integ.t), dt=$dt" minlen extrema(dϕ[B⁻])#, next at t=$(integ.t+dt)" 
#     if verbose
#         # reinit_err = sdf_err_L∞(ϕ, dom)
#         # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
#         @info "Reinit at t=$(integ.t), dt=$dt" dryfrac Tf Tw integ.t/dryfrac #, next at t=$(integ.t+dt)" 
#     end

#     return integ.t + dt
# end

"""
    sim_heatonly(fullconfig; tf=1e5, verbose=false)

Given a simulation configuration `fullconfig`, run a simulation.

This simulation is stripped-down: no mass transfer, no variation in ice & glass temperature

Maximum simulation time is specified by `tf`.
`verbose=true` will put out some info messages about simulation progress, i.e. at each reinitialization.

`fullconfig` should have the following fields:
- `init_prof`, types listed for [`make_ϕ0`](@ref)
- `fillvol`, fill volume with Unitful units (e.g. `3u"mL"`)
- `vialsize`, a string (e.g. `"2R"`) giving vial size
- `simgridsize`, a tuple/arraylike giving number of grid points to use for simulation. Defaults to `(51, 51)`.
- `Tf0`, an ice temperature with Unitful units 
- `Tw0`, a glass temperature (if the same as Tf0, can leave this out)
- `controls`, which has following fields (either scalar or array, with same length as `t_samp`:
    - `t_samp`, sampled measurement times. Needed only if other measurements are given during time
    - `Tsh`, shelf temperature: either a scalar (constant for full time span) or an array at specified time, in which case implemented via callback
    - `Q_ic`, ice RF heating. 
- `cparams`, which in turn has fields with Unitful units
    - `Kw`, 
    - `Kv` : heat transfer coefficients shelf
    - `Q_ck` : volumetric heating in cake 
    - `k`: thermal conductivity of cake
    - `ρf`: density of ice
    - `ΔH` : heat of sublimation of ice
    - `ϵ` : porosity of porous medium

Other parameters will be ignored.

During simulation, at each value of `t_samp`, the values of any `controls` which are arrays will be added to an internal dict called `params`.

If you pass in an array of values for multiple of `Tsh`, `Q_gl_RF`, or others, they must all have the same length as `t_samp`.

If you are getting a warning about instability, it can sometimes be fixed by tinkering with the reinitialization behavior.


"""
function sim_heatonly(fullconfig; tf=1e5, verbose=false)

    # ------------------- Get simulation parameters

    @unpack cparams, init_prof, Tf0, controls, vialsize, fillvol = fullconfig

    # Default values for non-essential parameters
    Tw0 = get(fullconfig, :Tw0, Tf0) # Default to same ice & glass temperature if glass initial not given
    simgridsize = get(fullconfig, :simgridsize, (51,51))

    # --------- Set up simulation domain
    r_vial = get_vial_radii(vialsize)[1]
    z_fill = fillvol / π / r_vial^2

    rmax = ustrip(u"m", r_vial)
    zmax = ustrip(u"m", z_fill)

    dom = Domain(simgridsize..., rmax, zmax)

    # ----- Nondimensionalize everything

    Tf0 = ustrip(u"K", Tf0)
    Tw0 = ustrip(u"K", Tw0)
    params, ncontrols = params_nondim_setup(cparams, controls) # Covers the various physical parameters 
    # if verbose
    #     @info "Variables used in callback:" meas_keys
    # end
    


    ϕ0 = make_ϕ0(init_prof, dom)
    if verbose
        @info "Initializing ϕ"
    end
    # Make sure that the starting profile is very well-initialized
    # The chosen tolerance is designed to the error almost always seen in norm of the gradient
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=10000, tol=1.2/max(dom.nr,dom.nz), err_reg=:all) 

    ϕ0_flat = reshape(ϕ0, :)

    
    # Full array of starting state variables ------------
    # u0 = similar(ϕ0_flat, dom.ntot+2) # Add 2 to length: Tf, Tw
    # u0[1:dom.ntot] .= ϕ0_flat
    # u0[dom.ntot+1] = Tf0 
    # u0[dom.ntot+2] = Tw0 
    u0 = make_u0_ndim(init_prof, Tf0, Tw0, dom)

    # Cached array for using last pressure state as guess
    p_last = fill(0.0, size(dom))

    # ----- Set up parameters dictionary and measurement callback
    meas_affect!(integ) = input_measurements!(integ, meas_keys, ncontrols)
    cb_meas = PresetTimeCallback(ncontrols[:t_samp], meas_affect!, filter_tstops=true)

    # ---- Set up ODEProblem
    prob_pars = (dom, params, p_last)
    tspan = (0, tf)
    prob = ODEProblem(dudt_heatonly!, u0, tspan, prob_pars)

    # --- Set up reinitialization callback

    cb_reinit = IterativeCallback(x->next_reinit_time_heatonly(x, verbose=verbose), reinit_wrap,  initial_affect = true)

    # --- Set up simulation end callback

    # When the minimum value of ϕ is 0, front has disappeared
    cond_end(u, t, integ) = minimum(u) 
    # ContinuousCallback gets thrown when `cond` evaluates to 0
    # `terminate!` ends the solve there
    cb_end = ContinuousCallback(cond_end, terminate!)

    # ------- Put together callbacks 
    cbs = CallbackSet(cb_reinit, cb_end, cb_meas)

    if verbose
        @info "Beginning solve"
    end
    # --- Solve
    sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping: default
    # sol = solve(prob, SSPRK33(), dt=1e-4, callback=cbs; ) # Fixed timestepping
    # sol = solve(prob, Tsit5(), callback=cbs; ) # Different adaptive integrator
    return @strdict sol dom
end