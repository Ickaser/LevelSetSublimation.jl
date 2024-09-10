# -------------------------- Heat transfer only functions
# export compute_frontvel_heat, dudt_heatonly!
# export compute_Qice


"""
    dudt_heatonly!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

This function leaves `Tf` and `Tw` untouched, since there isn't a way to govern their dynamics without mass transfer.

`u` and `du` are both structured as follows:
First `dom.ntot` values are `ϕ`, reshaped; `dom.ntot+1` index is frozen temperature `Tf`, `dom.ntot+2` index is glass temperature `Tw`

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function dudt_heatonly!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)
    ϕ, Tf, Tw = ϕ_T_from_u(u, dom)

    T = solve_T(u, Tf, dom, params)

    vf, dϕdx_all = compute_frontvel_heat(u, Tf, T, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    dTf .= 0
    dTw .= 0

    # dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)
        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind]>0, vz[ind]>0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end
    # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    # @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) Tf[1] params[:Tsh] Tw extrema(vr) extrema(vz)
    return nothing
end


"""
    dudt_heatonly(u, dom::Domain, params)
    dudt_heatonly(u, config)
    
Compute the time derivative of `u` with given parameters.

Wraps a call on `dudt_heatonly!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatonly(u, dom::Domain, params)
    integ_pars = (dom, params, zeros(Float64, size(dom)))
    du = similar(u)
    du[dom.ntot+1:end] .= 0
    dudt_heatonly!(du, u, integ_pars, 0.0)
    return du
end
function dudt_heatonly(u, config)
    @unpack vialsize, fillvol = config

    dom = Domain(config)
    params, ncontrols = params_nondim_setup(config[:cparams], config[:controls])

    dudt_heatonly(u, dom, params)
end


"""
    compute_frontvel_heat(u, T, p, dom::Domain, params; debug=false)

Generate an empty velocity field and compute velocity on `Γ⁺` (i.e. cells on Γ with ϕ>0). 
"""
function compute_frontvel_heat(u, Tf, T, dom::Domain, params; debug=false)

    @unpack kd, ΔH, ρf, ϵ = params
    ϕ = ϕ_T_from_u(u, dom)[1]

    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
     
    vf = zeros(eltype(ϕ), dom.nr, dom.nz, 2)
    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    Qice = compute_Qice_nodry(u, T, dom, params)
    Q_ice_per_surf = Qice / compute_icesurf_δ(ϕ, dom)

    for c in Γ⁺
        ir, iz = Tuple(c)
        dTr, dTz = compute_Tderiv(u, Tf, T, ir, iz, dom, params)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, dTr<0, dTz<0, dϕdx_all, dom)

        # Normal is out of the ice
        # md >0 for sublimation occurring
        # v = md/ρ * -∇ϕ
        q_grad = kd* (dTr*dϕdr + dTz * dϕdz) 
        md_l = (q_grad + Q_ice_per_surf)/ΔH
        vtot = md_l / ρf / ϵ 
        vf[c,1] = -vtot * dϕdr
        vf[c,2] = -vtot * dϕdz

    end
    
    return vf, dϕdx_all
end

"""
    function compute_Qice(ϕ, T, p, dom::Domain, params)

Compute the total heat input into frozen & dried domains. Also passes glass-ice heat as separate return.

See p. 112-113 from paperlike notes. At pseudosteady conditions, all heat getting added to the system goes to the frozen domain,
so we don't actually need to treat different areas of the surface.
In addition, the sublimation flux is simply evaluated at the top cake surface.

"""
function compute_Qice(u, T, p, dom::Domain, params)

    @unpack ΔH= params
    # ϕ, Tf, Tw = ϕ_T_from_u(u, dom)[1]

    Q_glshvol, Q_vwf = compute_Qice_noflow(u, T, dom, params)

    # Sublimation rate
    md = compute_topmassflux(u, T, p, dom, params)
    Qsub = - md * ΔH

    return Q_glshvol + Qsub, Q_vwf
end


"""
    compute_Qice_noflow(u, T, dom::Domain, params)

Compute the total heat input into frozen & dried domains, assuming mass flow is zero. Also passes glass-ice heat as separate return.
"""
function compute_Qice_noflow(u, T, dom::Domain, params)

    @unpack Kshf, Kvwf, QRFf, Q_ck, Tsh= params
    ϕ, Tf, Tw = ϕ_T_from_u(u, dom)

    # Heat flux from shelf, at bottom of vial
    rweights = ones(Float64, dom.nr) 
    rweights[begin] = rweights[end] = 0.5
    Qsh = 2π* Kshf * dom.dr* sum(rweights .* dom.rgrid .* (Tsh .- T[:,1] ))

    # Heat flux from glass, at outer radius
    zweights = ones(Float64, dom.nz) 
    zweights[begin] = zweights[end] = 0.5
    Q_vwf = 2π*dom.rmax * Kvwf * dom.dz * sum(zweights .* ( Tw .- T[end,:]))

    # Volumetric heat in cake and ice
    icevol = compute_icevol(ϕ, dom)
    dryvol = π*dom.rmax^2*dom.zmax - icevol
    Qvol = icevol * QRFf + dryvol * Q_ck
    
    # @info "Q" Tsh Tf Qsh Q_vwf Qvol icevol

    return Qsh + Q_vwf + Qvol, Q_vwf
end

"""
    compute_Qice_nodry(u, T, dom::Domain, params)

Compute the total heat input into frozen domain from volumetric, shelf, and glass. 
Contrast with `compute_Qice_noflow` and `compute_Qice`, which include heat to dried domain.
"""
function compute_Qice_nodry(u, T, dom::Domain, params)
    @unpack Kshf, Kvwf, QRFf, Q_ck, Tsh = params
    ϕ, Tf, Tw = ϕ_T_from_u(u, dom)

    # Heat flux from shelf, at bottom of vial
    rweights = compute_icesh_area_weights(ϕ, dom)
    Qsh = 2π* Kshf * sum(rweights .* (Tsh .- T[:,1] ))
    

    # Heat flux from glass, at outer radius
    zweights = compute_icegl_area_weights(ϕ, dom)
    Q_vwf = 2π*dom.rmax * Kvwf * sum(zweights .* ( Tw .- T[end,:]))

    # Volumetric heat in cake and ice
    icevol = compute_icevol(ϕ, dom)
    Qvol = icevol * QRFf 

    # @info "geom" icevol get_subf_r(ϕ, dom)
    return Qsh + Q_vwf + Qvol
end

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
    - `QRFf`, ice RF heating. 
- `cparams`, which in turn has fields with Unitful units
    - `Kvwf`, 
    - `Kshf` : heat transfer coefficients shelf
    - `Q_ck` : volumetric heating in cake 
    - `kd`: thermal conductivity of cake
    - `ρf`: density of ice
    - `ΔH` : heat of sublimation of ice
    - `ϵ` : porosity of porous medium

Other parameters will be ignored.

During simulation, at each value of `t_samp`, the values of any `controls` which are arrays will be added to an internal dict called `params`.

If you pass in an array of values for multiple of `Tsh`, `QRFvw`, or others, they must all have the same length as `t_samp`.

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