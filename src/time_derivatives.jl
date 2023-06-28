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
    Tf_last = integ_pars[4]
    dϕ, dTf, dTgl = ϕ_T_from_u_view(du, dom)
    ϕ, Tf, Tgl = ϕ_T_from_u_view(u, dom)
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    # if any(Tf .!= clamp.(Tf, 200, 350))
    #     @warn "Crazy Tf, clamped" Tf
    #     clamp!(Tf, 200, 350)
    # end
    # if Tgl != clamp(Tgl, 200, 400)
    #     @warn "Crazy Tgl, clamped"
    #     Tgl = clamp(Tgl, 200, 400)
    # end
    # T = solve_T(u, dom, params)

    # p_sub = calc_psub.(Tf)
    # # If no sublimation occurring, just let temperature increase
    # if all(p_sub .< params[:p_ch]) # No driving force for mass transfer: no ice loss, just temperature change
    #     Qice, Qgl = compute_Qice_noflow(u, T, dom, params)
    #     # TODO 
    #     dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-8) # Prevent explosion during last time step by not letting volume go to 0
    #     dTgl .= (Q_gl_RF - Qgl) / m_cp_gl
    #     dϕ .= 0.0
    #     return nothing
    # end

    # p = solve_p(u, T, dom, params, p_last)

    Tfs, T, p = pseudosteady_Tf_T_p(u, dom, params, Tf_last, p_last)

    # Drive guess in the direction of the solved value, with approximate dt
    dTf .= (Tfs - Tf) /60 
    # Store pseudosteady in u, for use in other functions
    Tf .= Tfs

    integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    integ_pars[4] .= Tfs
    vf, dϕdx_all = compute_frontvel_mass(u, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # TODO
    Qice, Qgl = compute_Qice(u, T, p, dom, params)
    # dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    # dTfdt_radial!(dTf, u, T, p, dϕdx_all, dom, params)
    dTgl .= (Q_gl_RF - Qgl) / m_cp_gl

    # dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)

        # # Boundary cases: use internal derivative
        # if ir == dom.nr # Right boundary
        #     dϕdr = dϕdr_w[ind]
        # elseif ir == 1 # Left boundary
        #     dϕdr = dϕdr_e[ind]
        # else
        #     dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
        # end
        # # Boundary cases: use internal derivative
        # if iz == dom.nz # Top boundary
        #     dϕdz = dϕdz_s[ind]
        # elseif iz == 1 # Bottom boundary
        #     dϕdz = dϕdz_n[ind]
        # else
        #     dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])
        # end

        # # Don't treat boundaries differently
        # dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
        # dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end
    dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) Tgl
    if minimum(dϕ) < 0
        @info "negative dϕ" findall(dϕ.<0)
    end
    return nothing
end

function choose_dϕdx_boundary(ir, iz, west_true::Bool, south_true::Bool, dϕdx_all, dom::Domain)
    # Boundary cases: use internal derivative
    if ir == dom.nr # Right boundary
        dϕdr = dϕdx_all[1][ir, iz]
    elseif ir == 1 # Left boundary
        dϕdr = dϕdx_all[2][ir,iz]
    else
        dϕdr = (west_true > 0 ? dϕdx_all[1][ir,iz] : dϕdx_all[2][ir,iz])
    end
    # Boundary cases: use internal derivative
    if iz == dom.nz # Top boundary
        dϕdz = dϕdx_all[3][ir, iz]
    elseif iz == 1 # Bottom boundary
        dϕdz = dϕdx_all[4][ir, iz]
    else
        dϕdz = (south_true ? dϕdx_all[3][ir, iz] : dϕdx_all[4][ir, iz])
    end

    # # Don't treat boundaries differently
    # dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
    # dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])
    return dϕdr, dϕdz
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

# ---------------------------
function local_sub_heating_dϕdx(u, T, p, ir, iz, dϕdx_all, dom, params)
    @unpack k, ΔH = params
    b = eval_b_loc(T, p, ir, iz, params)

    dTdr, dTdz = compute_Tderiv(u, T, ir, iz, dom, params)
    dpdr, dpdz = compute_pderiv(u, T, p, ir, iz, dom, params)
    # dϕdr = dpdr < 0 ? dϕdx_all[1][ir, iz] : dϕdx_all[2][ir, iz] # East or west based on mass flow
    # dϕdz = dpdz < 0 ? dϕdx_all[3][ir, iz] : dϕdx_all[4][ir, iz] # North or south based on mass flow
    dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, dpdr<0, dpdz<0, dϕdx_all, dom)

    qflux = k*(dϕdr*dTdr + dϕdz*dTdz)
    mflux = b*(dϕdr*dpdr + dϕdz*dpdz)
    return qflux + ΔH*mflux, dϕdr, dϕdz
end

function dTfdt_radial(u, T, p, dϕdx_all, dom::Domain, params)
    dTfdt = similar(u, dom.nr)
    dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom, params)
    return dTfdt
end

function dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom::Domain, params)
    ϕ, Tgl = ϕ_T_from_u(u, dom)[[true, false, true]]
    @unpack ρf, Cpf, kf, Q_ic = params

    Δξ, bot_contact, top_contact = compute_iceht_bottopcont(ϕ, dom)

    has_ice = (Δξ .> 0)
    no_ice = (Δξ .== 0)

    for ir in axes(dTfdt, 1)
        if no_ice[ir]
            continue
        end

        if ir == 1
            dTfdr = 0
            d2Tfdr2 = (-2Tf[1] + 2Tf[2])*dom.dr2 # Adiabatic ghost cell
        elseif ir == dom.nr
            dTfdr = params[:Kgl]/kf*(Tgl - Tf[ir])
            d2Tfdr2 = (-2Tf[dom.nr] + 2Tf[dom.nr-1] + 2*dom.dr*dTfdr)*dom.dr2 # Robin ghost cell
        elseif no_ice[ir-1] # On left side: away from center
            @error "Not implemented"
        elseif no_ice[ir+1] # On right side: pulled away from wall
            integ_cells = [CI(ir+1, iz) for iz in 1:dom.nz if ϕ[ir,iz]<=0]
            surf_integral = 0
            for cell in integ_cells
                q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, T, p, Tuple(cell)..., dϕdx_all, dom, params)
                δ = compute_local_δ(cell, ϕ, dom )
                surf_integral += dom.rgrid[ir+1]*q*δ*dom.dr*dom.dz / dϕdr
            end

            dTfdr = 1/kf/dom.rgrid[ir]/Δξ[ir] * surf_integral
            # Use a ghost cell
            d2Tfdr2 = (-2Tf[ir] + 2Tf[ir-1] + 2*dom.dr*dTfdr)*dom.dr2
        else
            dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        end

        if top_contact[ir]
            # Adiabatic boundary
            top_bound_term = 0
        else
            # Stefan boundary
            iz = findlast(ϕ[ir,:] .<=0) + 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, T, p, ir, iz, dϕdx_all, dom, params)
            top_bound_term = q - kf*dϕdr/dϕdz*dTfdr
        end

        if bot_contact[ir]
            # Shelf boundary
            bot_bound_term = params[:Kv]*(T[ir,begin] - params[:Tsh])
        else
            # Stefan boundary
            iz = findfirst(ϕ[ir,:] .<=0) - 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, T, p, ir, iz, dϕdx_all, dom, params)
            bot_bound_term = q - kf*dϕdr/dϕdz*dTfdr
        end

        # @info "vals" dTfdr d2Tfdr2 Q_ic top_bound_term bot_bound_term 
        dTfdt[ir] = (kf*((ir == 1 ? 0 : 1/dom.rgrid[ir])*dTfdr + d2Tfdr2) + Q_ic + 
            (top_bound_term - bot_bound_term)/Δξ[ir])/ρf/Cpf
    end

    dTfdt[no_ice] .= 0 # Set to 0 elsewhere

    # TODO
    # Need to define Tf one cell beyond the ice, actually, for other parts of implementation.
    # So: set up a ghost cell, with dTfdr defined by Stefan boundary
    # Question is: how to evaluate, exactly? Average across vertical direction?

    # Perhaps can use a fictitious dTfdt at that point to define the ghost cell

    lbound = [i for i in 2:dom.nr if (has_ice[i] && no_ice[i-1])]
    rbound = [i for i in 1:dom.nr-1 if (has_ice[i] && no_ice[i+1])]
    if length(lbound) == 0
        #nothing
    elseif length(lbound) == 1
        # Treat ghost cell
        # dTfdt[lbound[1]-1] = ...
    else
        @warn "Multiple chunks of ice may not be handled correctly." lbound rbound
    end
    if length(rbound) == 0
        #nothing
    elseif length(rbound) == 1
        # Treat ghost cell
        # dTfdt[rbound[1]+1] = ...
    else
        @warn "Multiple chunks of ice may not be handled correctly." lbound rbound
    end


    # Alternatively: freeze Tf values that no longer have ice, use them as is.
    # That seems like a bad idea though.



    # wallBC = params[:Kgl]/kf*(Tgl - Tf[ir])

    # dTfdr = fill(0.0, dom.nr)
    # dTfdr[dom.nr] = no_ice[dom.nr] ? 0 : wallBC
    # for ir in 2:dom.nr-1
    #     dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
    # end
    # d2Tfdr2 = fill(0.0, dom.nr)
    # # dTfdr[1] = 0
    # # dTfdr[dom.nr] = no_ice[dom.nr] ? 0 : params[:Kgl]/kf*(Tgl - Tf[ir])
    # d2Tfdr2[1] = (-2Tf[1] + 2Tf[2])*dom.dr2 # Adiabatic ghost cell
    # d2Tfdr2[dom.nr] = (-2Tf[dom.nr] + 2Tf[dom.nr-1] + 2*dom.dr*wallBC)*dom.dr2 # Adiabatic ghost cell
    # for ir in 2:dom.nr-1
    #     d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
    # end


    # top_bound_term = fill(0.0, dom.nr)
    # bot_bound_term = fill(0.0, dom.nr)
    # for ir in axes(ϕ, 1)
    #     if no_ice[ir]
    #         continue
    #     end
    #     if top_contact[ir]
    #         # Adiabatic boundary
    #         top_bound_term = 0
    #     else
    #         # Stefan boundary
    #         q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, T, p, ir, iz, dϕdx_all, dom, params)
    #         top_bound_term = q - kf*dϕdr/dϕdz*dTfdr[ir]
    #     end
    #     if bot_contact[ir]
    #         # Shelf boundary
    #         bot_bound_term = params[:Kv]*(T[ir,begin] - params[:Tsh])
    #     else
    #         # Stefan boundary
    #         q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, T, p, ir, iz, dϕdx_all, dom, params)
    #         bot_bound_term = q - kf*dϕdr/dϕdz*dTfdr[ir]
    #     end
    # end



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

    # dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)
        # Boundary cases: use internal derivative
        # if ir == dom.nr # Right boundary
        #     dϕdr = dϕdr_w[ind]
        # elseif ir == 1 # Left boundary
        #     dϕdr = dϕdr_e[ind]
        # else
        #     dϕdr = (vr[ind] > 0 ? dϕdr_w[ind] : dϕdr_e[ind])
        # end
        # # Boundary cases: use internal derivative
        # if iz == dom.nz # Top boundary
        #     dϕdz = dϕdz_s[ind]
        # elseif iz == 1 # Bottom boundary
        #     dϕdz = dϕdz_n[ind]
        # else
        #     dϕdz = (vz[ind] > 0 ? dϕdz_s[ind] : dϕdz_n[ind])
        # end
        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind]>0, vz[ind]>0, dϕdx_all, dom)

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