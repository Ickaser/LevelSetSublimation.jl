export dudt_heatmass, dudt_heatonly

"""
    dudt_heatmass!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

Splitting `u` and `du` into `ϕ`, `Tf`, and `Tgl` is handled by `ϕ_T_from_u` and `ϕ_T_from_u_view`.

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function dudt_heatmass!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    p_last = integ_pars[3]
    dϕ, dTf, dTgl = ϕ_T_from_u_view(du, dom)
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    if any(Tf .!= clamp.(Tf, 200, 350))
        @warn "Crazy Tf, clamped"
        clamp!(Tf, 200, 350)
    end
    if Tgl != clamp(Tgl, 200, 400)
        @warn "Crazy Tgl, clamped"
        Tgl = clamp(Tgl, 200, 400)
    end
    T = solve_T(u, dom, params)

    p_sub = calc_psub.(Tf)
    # If no sublimation occurring, just let temperature increase
    if all(p_sub .< params[:p_ch]) # No driving force for mass transfer: no ice loss, just temperature change
        Qice, Qgl = compute_Qice_noflow(u, T, dom, params)
        # TODO 
        dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-8) # Prevent explosion during last time step by not letting volume go to 0
        dTgl .= (Q_gl_RF - Qgl) / m_cp_gl
        dϕ .= 0.0
        return nothing
    end

    p = solve_p(u, T, dom, params, p_last)
    integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    vf, dϕdx_all = compute_frontvel_mass(u, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # TODO
    Qice, Qgl = compute_Qice(u, T, p, dom, params)
    dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    dTgl .= (Q_gl_RF - Qgl) / m_cp_gl
    # du[ntot+1] .= dTfdt
    # du[ntot+2] .= dTgldt

    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)

        # Boundary cases: use internal derivative
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

        # # Don't treat boundaries differently
        # dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
        # dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end
    # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    # @info "prog: t=$t, dryfrac=$dryfrac"
    return nothing
end


"""
    dudt_heatmass(u, dom::Domain, params)
    dudt_heatmass(u, config)
    
Compute the time derivative of `u` with given parameters.

`u` has `dom.ntot` entries for `ϕ`, `dom.nr` for `Tf`, and 1 for `Tgl`.

Wraps a call on `dudt_heatmass!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatmass(u, dom::Domain, params)
    integ_pars = (dom, params, zeros(Float64, size(dom)))
    du = similar(u)
    dudt_heatmass!(du, u, integ_pars, 0.0)
    return du
end
function dudt_heatmass(u, config)
    # Set up simulation domain & parameters
    dom = Domain(config)
    params, ncontrols = params_nondim_setup(config[:cparams], config[:controls])

    dudt_heatmass(u, dom, params)
end

"""
    dudt_heatmass_params(u, config)
    
Compute the time derivative of `u` with given parameters, and also return `dom` and `params` associated with the given `config`.

Wraps a call on `dudt_heatmass`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatmass_params(u, config)
    # Set up simulation domain & parameters
    dom = Domain(config)
    params, ncontrols = params_nondim_setup(config[:cparams], config[:controls])

    dudt_heatmass(u, dom, params), dom, params
end

# -------------------------- Heat transfer only functions


"""
    dudt_heatonly!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

This function leaves `Tf` and `Tgl` untouched, since there isn't a way to govern their dynamics without mass transfer.

`u` and `du` are both structured as follows:
First `dom.ntot` values are `ϕ`, reshaped; `dom.ntot+1` index is frozen temperature `Tf`, `dom.ntot+2` index is glass temperature `Tgl`

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function dudt_heatonly!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    p_last = integ_pars[3]
    dϕ, dTf, dTgl = ϕ_T_from_u_view(du, dom)
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)

    T = solve_T(u, dom, params)

    vf, dϕdx_all = compute_frontvel_heat(u, T, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    dTf .= 0
    dTgl .= 0

    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)
        # Boundary cases: use internal derivative
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

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end
    # dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    # @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) Tf[1] params[:Tsh] Tgl extrema(vr) extrema(vz)
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