export sim_from_dict

export ϕevol_RHS, ϕ_T_from_u

# --------- Convenience functions that need a home

"""
    ϕ_T_from_u(u, dom)

Take the current system state `u` and break it into `ϕ`, `Tf`, and `Tgl`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u(u, dom)
    ϕ = reshape(u[1:dom.ntot], size(dom))
    Tf = u[dom.ntot+1]
    Tgl = u[dom.ntot+2]
    return ϕ, Tf, Tgl
end

# ---------- Fully adaptive time stepping functions

"""
    ϕevol_RHS!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

`u` and `du` are both structured as follows:
First `dom.ntot` values are `ϕ`, reshaped; `dom.ntot+1` index is frozen temperature `Tf`

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function ϕevol_RHS!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    p_last = integ_pars[3]
    ntot = dom.ntot
    dϕ = reshape((@view du[1:ntot]), dom.nr, dom.nz)
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    u[dom.ntot+1] = clamp(Tf, 200, 350)  # Prevent crazy temperatures from getting passed through to other functions
    u[dom.ntot+2] = clamp(Tgl, 200, 400) # Prevent crazy temperatures from getting passed through to other functions
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    T = solve_T(u, dom, params)

    p_sub = calc_psub(Tf) 
    if p_sub < params[:p_ch] # No driving force for mass transfer: no ice loss, just temperature change
        Qice, Qgl = compute_Qice_noflow(u, T, dom, params)
        dTfdt = Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
        dTgldt =  (Q_gl_RF - Qgl) / m_cp_gl
        dϕ .= 0.0
        du[ntot+1] = dTfdt
        du[ntot+2] = dTgldt
        # @info "prog nosub"  Tgl Qgl Q_gl_RF dTgldt t
        return nothing
    end

    p = solve_p(u, T, dom, params, p0 = p_last)
    integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    vf, dϕdx_all = extrap_v_fastmarch(u, T, p, dom, params)
    vr = @view vf[:,:,1]
    vz = @view vf[:,:,2]
    # pl1 = heat(p, dom)
    # plot_contour(ϕ, dom)
    # pl2 = heat(vz, dom)
    # plot_contour(ϕ, dom)
    # display(plot(pl1, pl2))
    # @info "test" p[25:35,40:51] vr[25:35,40:51] vz[25:35,40:51]
    
    Qice, Qgl = compute_Qice(u, T, p, dom, params)
    dTfdt = Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    dTgldt = (Q_gl_RF - Qgl) / m_cp_gl
    du[ntot+1] = dTfdt
    du[ntot+2] = dTgldt
    # @info "glass heating?" dTgldt Qgl

    # Compute dϕ/dt = - v ⋅ ∇ ϕ
    # indmin = CI(1, 1)
    # indmax = CI(dom.nr, dom.nz)
    # rshift = [CI(i, 0) for i in -3:3]
    # zshift = [CI(0, i) for i in -3:3]
    # dϕdr_e = zeros(Float64, size(dom))
    # dϕdr_w = zeros(Float64, size(dom))
    # dϕdz_n = zeros(Float64, size(dom))
    # dϕdz_s = zeros(Float64, size(dom))
    dϕdr_e, dϕdr_w, dϕdz_n, dϕdz_s = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)
    #     rst = max.(min.([ind].+rshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
    #     zst = max.(min.([ind].+zshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
    #     dϕdr_w[ind], dϕdr_e[ind] = wenodiffs_local(ϕ[rst]..., dom.dr)
    #     dϕdz_s[ind], dϕdz_n[ind] = wenodiffs_local(ϕ[zst]..., dom.dz)
    #     # Boundary cases: use internal derivative
        if ir == dom.nr # Right boundary
            dϕdr = dϕdr_w[ind]
        elseif ir == 1 # Left boundary
            dϕdr = dϕdr_e[ind]
        else
            dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
        end
        # Boundary cases: use internal derivative
        if iz == dom.nz # Top boundary
            dϕdz = dϕdz_s[ind]
        elseif iz == 1 # Bottom boundary
            dϕdz = dϕdz_n[ind]
        else
            dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])
        end

        rcomp = dϕdr*vr[ind]
        zcomp = dϕdz*vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp 
    end
    dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    # @info "prog: t=$t, dryfrac=$dryfrac" maximum(dϕ)
    return nothing
end


"""
    ϕevol_RHS(u, dom::Domain, params)
    ϕevol_RHS(u, config)
    
Compute the time derivative of `u` with given parameters.

`u` has `dom.ntot` entries for `ϕ` and one each for `Tf` and `Tgl`.

Wraps a call on `ϕevol_RHS!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function ϕevol_RHS(u, dom::Domain, params)
    integ_pars = (dom, params, zeros(Float64, size(dom)))
    du = similar(u)
    # dϕ = zeros(dom.nr, dom.nz)

    # dϕ_flat = reshape(dϕ, :)
    # u = similar(ϕ, dom.ntot+1)
    # u[1:dom.ntot] .= reshape(ϕ, :)
    # u[dom.ntot+1] = params[:Tf]
    ϕevol_RHS!(du, u, integ_pars, 0.0)
    return du
end
function ϕevol_RHS(u, config, t=0)
    @unpack vialsize, fillvol = config
    simgridsize = get(config, :simgridsize, (51,51))

    # --------- Set up simulation domain
    r_vial = get_vial_rad(vialsize)
    z_fill = fillvol / π / r_vial^2

    rmax = ustrip(u"m", r_vial)
    zmax = ustrip(u"m", z_fill)

    dom = Domain(simgridsize..., rmax, zmax)

    params, ncontrols = params_nondim_setup(config[:cparams], config[:controls])

    ϕevol_RHS(u, dom, params)
end

"""
    reinit_wrap(integ)

Thin wrapper to reinitialize the state of the level set function.

Calls `reinitialize_ϕ!(ϕ, dom)`, so uses the default reinitialization setup.
Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
"""
function reinit_wrap(integ; verbose=false)
    if verbose
        @info "Reinit at t=$(integ.t)"
    end
dom = integ.p[1]
ϕ = reshape((@view integ.u[1:dom.ntot]), dom.nr, dom.nz)
# reinitialize_ϕ!(ϕ, dom) 
ϕa = copy(ϕ)
reinitialize_ϕ_HCR!(ϕ, dom, maxsteps=100, tol=1/max(dom.nr, dom.nz)) 
# pl1 = Plots.contour(dom.rgrid, dom.zgrid, ϕa', aspect_ratio=:equal)
# pl2 = Plots.contour(dom.rgrid, dom.zgrid, ϕ' , aspect_ratio=:equal)
# pl3 = heat(𝒢_weno_all(ϕa, dom), dom)
# plot_contour(ϕa, dom)
# pl4 = heat(𝒢_weno_all(ϕ, dom), dom)
# plot_contour(ϕ, dom)
# # display(plot(pl1, pl2))
# savefig(plot(pl1, pl2, pl3, pl4, size=(800,800)), plotsdir("reinit_testing", "reinit_snap_$(round(integ.t)).png"))
end

"""
next_reinit_time(integ)

Compute the next time reinitialization should be necessary, given the current integrator state.
Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
"""
function next_reinit_time(integ; verbose=false)
    dom = integ.p[1]
    du = similar(integ.u)
    ϕevol_RHS!(du, integ.u, integ.p, integ.t)

    # The main region of concern is the frozen region near interface
    # Find the largest value of dϕdt in that region
    ϕ, Tf, Tgl = ϕ_T_from_u(integ.u, dom)
    dϕ, dTfdt, dTgldt = ϕ_T_from_u(du, dom)
    B = identify_B(ϕ, dom)
    B⁻ = B .& (ϕ .<= 0)
    # @info "something's up" ϕ dϕ B sum(B⁻)
    if sum(B⁻) > 0
        max_dϕdt = maximum(abs.(dϕ[B⁻]))
    else
        max_dϕdt = maximum(abs.(dϕ))
    end
    @unpack p_ch = integ.p[2]
    # if max_dϕdt == 0 # If no sublimation is happening...
    if calc_psub(Tf) < p_ch
        # Guess when it will start
        # Find temperature such that p_sub = p_ch
        Tstart = Tf
        for i in 1:5
            pg = calc_psub(Tstart)
            dpdT = pg - calc_psub(Tstart - 1)
            dT = -(pg-p_ch)/dpdT
            Tstart += dT
            # @info "Newton iter" Tstart i
        end
        dt = abs((Tstart - Tf)/dTfdt) * 1.05 # Overshoot slightly, to get into real drying behavior before next check time
        # @info "Sublimation stopped, estimating reinit." calc_psub(Tstart) calc_psub(Tf) p_ch dt dTfdt
        return integ.t + dt
    end

    # Reinit next when interface should have moved across half of band around interface 
    domfrac = 0.1 # Should reinit roughly 10-15 times during simulation
    minlen = min(domfrac*dom.rmax , domfrac*dom.zmax)  # Also: at early times, do more often
    dt = minlen / max_dϕdt
    dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    if dryfrac < domfrac || integ.t < 1000
        # dt = dt^(1//3) # Reinit more often in early stages
        dt = clamp(dt, 1, 100)
    end
    # dt = 0.5 * minlen / max_dϕdt 
    # @info "Reinit at t=$(integ.t), dt=$dt" minlen extrema(dϕ[B⁻])#, next at t=$(integ.t+dt)" 
    if verbose
        reinit_err = sdf_err_L∞(ϕ, dom)
        # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "Reinit at t=$(integ.t), dt=$dt" dryfrac Tf Tgl reinit_err #, next at t=$(integ.t+dt)" 
    end

    return integ.t + dt
end

function cond_reinit(u, t, integ)
    dom = integ.p[1]
    ϕ = ϕ_T_from_u(u, size(dom))[1]
    err = sdf_err_L1(ϕ, dom)
    tol = 1e-5 # Roughly dx^4
    @info "error from sdf" err-tol
    return err-tol
end


function input_measurements!(integ, meas_keys::Vector, controls)
    t_samp = controls[:t_samp]
    ti = argmin(abs.(t_samp .- integ.t))
    if !(integ.t ≈ t_samp[ti])
        @error "Issue with time sampling" integ.t t_samp[ti]
    end
    params = integ.p[2]
    for key in meas_keys
        params[key] = controls[key][ti]
    end
end

function input_measurements!(integ, meas_keys::Nothing, controls)
end

"""
    sim_from_dict(fullconfig; tf=100, verbose=false)

Given a simulation configuration `fullconfig`, run a simulation.

Maximum simulation time is specified by `tf`.
`verbose=true` will put out some info messages about simulation progress, i.e. at each reinitialization.

`fullconfig` should have the following fields:
- `ϕ0type`, types listed for [`make_ϕ0`](@ref)
- `fillvol`, fill volume with Unitful units (e.g. `3u"mL"`)
- `vialsize`, a string (e.g. `"2R"`) giving vial size
- `simgridsize`, a tuple/arraylike giving number of grid points to use for simulation. Defaults to `(51, 51)`.
- `Tf0`, an initial ice temperature with Unitful units 
- `Tgl0`, an initial glass temperature (if the same as Tf0, can leave this out)
- `controls`, which has following fields (either scalar or array, with same length as `t_samp`:
    - `t_samp`, sampled measurement times. Needed only if other measurements are given during time
    - `Tsh`, shelf temperature: either a scalar (constant for full time span) or an array at specified time, in which case implemented via callback
    - `Q_gl_RF`, glass RF heating. Scalar or array, like Tsh
    - `Q_ic`, ice RF heating. 
    - `p_ch` : pressure at top of cake
- `cparams`, which in turn has fields with Unitful units
    - `Kgl`, 
    - `Kv` : heat transfer coefficients shelf
    - `Q_ck` : volumetric heating in cake 
    - `k`: thermal conductivity of cake
    - `m_cp_gl` total thermal mass of ice, relevant to heating/cooling of glass wall
    - `ρf`: density of ice
    - `Cpf`: heat capacity of ice
    - `ΔH` : heat of sublimation of ice
    - `ϵ` : porosity of porous medium
    - `l` : dusty gas model constant: characteristic length for Knudsen diffusion
    - `κ` : dusty gas model constant: length^2 corresponding loosely to Darcy's Law permeability
    - `Mw`: molecular weight of species (water), with appropriate units
    - `μ` : dynamic viscosity of species (water), with appropriate units

During simulation, at each value of `t_samp`, the values of any `controls` which are arrays will be added to an internal dict called `params`.

If you pass in an array of values for multiple of `Tsh`, `Q_gl_RF`, or others, they must all have the same length as `t_samp`.

If you are getting a warning about instability, it can sometimes be fixed by tinkering with the reinitialization behavior.


"""
function sim_from_dict(fullconfig; tf=1e5, verbose=false)

    # ------------------- Get simulation parameters

    @unpack cparams, ϕ0type, Tf0, controls, vialsize, fillvol = fullconfig

    # Default values for non-essential parameters
    Tgl0 = get(fullconfig, :Tgl0, Tf0) # Default to same ice & glass temperature if glass initial not given
    simgridsize = get(fullconfig, :simgridsize, (51,51))

    # --------- Set up simulation domain
    r_vial = get_vial_rad(vialsize)
    z_fill = fillvol / π / r_vial^2

    rmax = ustrip(u"m", r_vial)
    zmax = ustrip(u"m", z_fill)

    dom = Domain(simgridsize..., rmax, zmax)

    # ----- Nondimensionalize everything

    Tf0 = ustrip(u"K", Tf0)
    Tgl0 = ustrip(u"K", Tgl0)
    params, meas_keys, ncontrols = params_nondim_setup(cparams, controls) # Covers the various physical parameters 
    


    ϕ0 = make_ϕ0(ϕ0type, dom)
    if verbose
        @info "Initializing ϕ"
    end
    # Make sure that the starting profile is very well-initialized
    # The chosen tolerance is designed to the error almost always seen in norm of the gradient
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=10000, tol=1.25/max(dom.nr,dom.nz)) 
    ϕ0_flat = reshape(ϕ0, :)

    
    # Full array of starting state variables ------------
    u0 = similar(ϕ0_flat, dom.ntot+2) # Add 2 to length: Tf, Tgl
    u0[1:dom.ntot] .= ϕ0_flat
    u0[dom.ntot+1] = Tf0 
    u0[dom.ntot+2] = Tgl0 
    p_sub = calc_psub(Tf0)
    # @info "Initial Tf:" Tf0 p_sub

    # Cached array for using last pressure state as guess
    p_last = fill(p_sub, size(dom))

    # ----- Set up parameters dictionary and measurement callback
    meas_affect!(integ) = input_measurements!(integ, meas_keys, ncontrols)
    cb_meas = PresetTimeCallback(ncontrols[:t_samp], meas_affect!, filter_tstops=true)

    # ---- Set up ODEProblem
    prob_pars = (dom, params, p_last)
    tspan = (0, tf)
    prob = ODEProblem(ϕevol_RHS!, u0, tspan, prob_pars)

    # --- Set up reinitialization callback

    # cb1 = PeriodicCallback(reinit_wrap, reinit_time, initial_affect=true)
    cb_reinit = IterativeCallback(x->next_reinit_time(x, verbose=verbose), reinit_wrap,  initial_affect = true)
    # if verbose
    #     cb_reinit = IterativeCallback(x->next_reinit_time(x, verbose), x->reinit_wrap(x, verbose=true),  initial_affect = true)
    # else
    #     cb_reinit = IterativeCallback(next_reinit_time, reinit_wrap,  initial_affect = true)
    # end


    # --- Set up simulation end callback

    # When the minimum value of ϕ is 0, front has disappeared
    cond_end(u, t, integ) = minimum(u) 
    # ContinuousCallback gets thrown when `cond` evaluates to 0
    # `terminate!` ends the solve there
    cb_end = ContinuousCallback(cond_end, terminate!)

    # ------- Put together callbacks 
    cbs = CallbackSet(cb_reinit, cb_end, cb_meas)

    # --- Solve
    sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping: default
    # sol = solve(prob, SSPRK33(), dt=1e-4, callback=cbs; ) # Fixed timestepping
    # sol = solve(prob, Tsit5(), callback=cbs; ) # Different adaptive integrator
    return @strdict sol dom
end



