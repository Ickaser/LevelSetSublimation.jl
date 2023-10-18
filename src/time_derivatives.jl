export dudt_heatmass, dudt_heatonly

"""
    dudt_heatmass!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

Splitting `u` and `du` into `ϕ`, `Tf`, and `Tw` is handled by `ϕ_T_from_u` and `ϕ_T_from_u_view`.

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function dudt_heatmass!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    p_last = integ_pars[3]
    Tf_last = integ_pars[4]
    controls = integ_pars[5]

    input_measurements!(params, t, controls)

    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)

    dTf .= 0 # Just in case some weird stuff got left there

    ϕ, Tw = ϕ_T_from_u_view(u, dom)[[true, false, true]]
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    if minimum(ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[:l])/length(params[:l])
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[:Mw]/params[:R]/minT)
        flux = (calc_psub(minT) - params[:p_ch])/dom.zmax*b_ave
        dϕ .= flux/params[:ρf]
        # dϕ .= 0
        # dTw .= 0
        @info "no ice" extrema(dϕ)
        return nothing
    end

    # p = solve_p(u, T, dom, params, p_last)

    # Tf = pseudosteady_Tf_T_p(u, dom, params, Tf_last, p_last)
    Tf = pseudosteady_Tf(u, dom, params, Tf_last)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)

    # Store pseudosteady in u, for use in other functions
    # Tf .= Tfs
    # dTf = (Tf - Tfg)/60

    integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    integ_pars[4] .= Tf
    vf, dϕdx_all = compute_frontvel_mass(u, Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # Compute time derivatives for 
    # dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    # dTfdt_radial!(dTf, u, T, p, dϕdx_all, dom, params)
    Qgl = compute_Qgl(u, T, dom, params)
    dTw .= (Q_gl_RF - Qgl) / m_cp_gl

    # dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end
    dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) extrema(T) Tw[1] params[:Tsh]
    if minimum(dϕ) < 0
        @info "negative dϕ" spy(ϕ .< 0) spy(dϕ .< 0) Tf[end]-Tf[1]
        # pl1 = heat(vr, dom)
        # pl2 = heat(vz, dom)
        # display(plot(pl1, pl2))
    end
    # if maximum(dϕ) > 1
    #     @info extrema(vz) extrema(vr) extrema(T) extrema(p) extrema(Tf)
    #     display(heat(T, dom))
    # end
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
    # return dϕdr, dϕdz
    normed = hypot(dϕdr, dϕdz)
    return dϕdr/normed, dϕdz/normed
end


"""
    dudt_heatmass(u, dom::Domain, params)
    dudt_heatmass(u, config)
    
Compute the time derivative of `u` with given parameters.

`u` has `dom.ntot` entries for `ϕ`, `dom.nr` for `Tf`, and 1 for `Tw`.

Wraps a call on `dudt_heatmass!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatmass(u, dom::Domain, params, ncontrols)
    integ_pars = (dom, params, fill(100.0, size(dom)), fill(233.15, dom.nr), ncontrols)
    du = similar(u)
    dudt_heatmass!(du, u, integ_pars, 0.0)
    return du
end
function dudt_heatmass(u, config, t=0)
    # Set up simulation domain & parameters
    dom = Domain(config)
    params, ncontrols = params_nondim_setup(config[:cparams], config[:controls])
    input_measurements!(params, t, ncontrols)
    dudt_heatmass(u, dom, params, ncontrols)
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

function set_params_at_t!(params, controls)
end


# ---------------------------
# For pseudosteady radial temperature

function local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
    @unpack k, ΔH = params
    b = eval_b_loc(T, p, ir, iz, params)

    dTdr, dTdz = compute_Tderiv(u, Tf, T, ir, iz, dom, params)
    dpdr, dpdz = compute_pderiv(u, Tf, T, p, ir, iz, dom, params)
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
    ϕ, Tw = ϕ_T_from_u(u, dom)[[true, false, true]]
    @unpack ρf, Cpf, kf, Q_ic, p_ch, ΔH, Rp0 = params

    Δξ, bot_contact, top_contact = compute_iceht_bottopcont(ϕ, dom)

    has_ice = (Δξ .> 0)
    no_ice = (Δξ .== 0)

    for ir in axes(dTfdt, 1)
        if no_ice[ir]
            continue
        end

        # Treat top and bottom
        if top_contact[ir]
            # # adiabatic boundary 
            # top_bound_term = 0
            # direct sublimation boundary
            top_bound_term = ΔH*(p_ch - calc_psub(Tf[ir]))/Rp0
            if typeof(top_bound_term) <: AbstractFloat
            end
        else
            iz = findlast(ϕ[ir,:] .<=0) + 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            # top_bound_term = q - kf*dϕdr/dϕdz*dTfdr
            top_bound_term = q*dϕdz
        end

        if bot_contact[ir]
            # Shelf boundary
            bot_bound_term = params[:Kv]*(T[ir,begin] - params[:Tsh])
        else
            # Stefan boundary
            iz = findfirst(ϕ[ir,:] .<=0) - 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            # bot_bound_term = q - kf*dϕdr/dϕdz*dTfdr
            bot_bound_term = q*dϕdz
        end

        # Treat radial derivatives

        if ir == 1
            dTfdr = 0
            d2Tfdr2 = (-2Tf[1] + 2Tf[2])*dom.dr2 # Adiabatic ghost cell
        elseif ir == dom.nr
            dTfdr = params[:Kw]/kf*(Tw - Tf[ir])
            d2Tfdr2 = (-2Tf[dom.nr] + 2Tf[dom.nr-1] + 2*dom.dr*dTfdr)*dom.dr2 # Robin ghost cell
        elseif no_ice[ir-1] && no_ice[ir+1] # Ice on both sides
            # @warn "Ice surrounded by gap: ignore radial gradients, treat only vertical." has_ice[ir-1:ir+1] 
            dTfdr = 0
            d2Tfdr2 = 0
            # dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            # d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        elseif no_ice[ir-1] # On left side: away from center
            # @warn "Not implemented carefully. Double check this math" ir has_ice[ir-1:ir+1] sum(has_ice) 
            integ_cells = [CI(ir-1, iz) for iz in 1:dom.nz if ϕ[ir,iz]<=0]
            surf_integral = 0.0
            surf_area = 0.0
            for cell in integ_cells
                q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(cell)..., dϕdx_all, dom, params)
                δm = compute_local_δ(cell, ϕ, dom )
                δp = compute_local_δ(cell + CI(1,0), ϕ, dom )
                loc_area = (δp*dom.rgrid[ir+1] + δm*dom.rgrid[ir])*2π*dom.dz*dom.dr
                surf_area += loc_area
                surf_integral += q*loc_area / dϕdr
            end

            # dTfdr = 1/kf/dom.rgrid[ir]/Δξ[ir] * surf_integral
            # dTfdr = surf_integral/surf_area /kf#*dom.rgrid[ir+1]/dom.rgrid[ir] 
            # flux = surf_integral #+ top_bound_term*top_contact[ir] - bot_bound_term*bot_contact[ir] 
            # dTfdr = flux/(dom.rgrid[ir]*2π*Δξ[ir])/kf # Divide integral by gridpoint surface area
            # dTfdr = surf_integral/(dom.rgrid[ir]*2π*Δξ[ir])/kf # Divide integral by gridpoint surface area
            dTfdr = surf_integral/surf_area/kf # Divide integral by gridpoint surface area
            # Use a ghost cell
            d2Tfdr2 = (-2Tf[ir] + 2Tf[ir-1] + 2*dom.dr*dTfdr)*dom.dr2
        elseif no_ice[ir+1] # On right side: pulled away from wall
            integ_cells = [CI(ir+1, iz) for iz in 1:dom.nz if ϕ[ir,iz]<=0]
            surf_integral = 0.0
            surf_area = 0.0
            for cell in integ_cells
                q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(cell)..., dϕdx_all, dom, params)
                δp = compute_local_δ(cell, ϕ, dom )
                δm = compute_local_δ(cell - CI(1,0), ϕ, dom )
                Hp = compute_local_H(cell, ϕ, dom)
                loc_area = (δp*dom.rgrid[ir+1] + δm*dom.rgrid[ir])*2π*dom.dz*dom.dr
                # loc_area = dom.rgrid[ir]*2π*dom.dz
                # loc_area = δ*dom.dz*dom.dr 
                # @info "check" cell halfcell loc_area
                surf_area += loc_area
                # surf_integral += q*loc_area/dϕdr #- Q_ic*Hp*dom.rgrid[ir]*2π*dom.dr
                surf_integral += q*loc_area/dϕdr #+ Q_ic*(Hp*dom.rgrid[ir]*dom.dr*2π*dom.dz)
                # TODO: should there be a dϕdr here?
            end

            # dTfdr = surf_integral/surf_area /kf#*dom.rgrid[ir+1]/dom.rgrid[ir] 
            # @info "checking" dom.rgrid[ir]*2π*Δξ[ir] surf_area
            # flux = surf_integral #+ top_bound_term*top_contact[ir] - bot_bound_term*bot_contact[ir] 
            # dTfdr = flux/(dom.rgrid[ir]*2π*Δξ[ir])/kf # Divide integral by gridpoint surface area
            dTfdr = surf_integral/surf_area/kf # Divide integral by gridpoint surface area
            # Use a ghost cell
            d2Tfdr2 = (-2Tf[ir] + 2Tf[ir-1] + 2*dom.dr*dTfdr)*dom.dr2
        else
            dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        end


        r1 = (ir == 1 ? 0 : 1/dom.rgrid[ir])
        dTfdt[ir] = (kf*(r1*dTfdr + d2Tfdr2) + Q_ic + 
            (top_bound_term - bot_bound_term)/Δξ[ir])/ρf/Cpf

        isnan(dTfdt[ir]) && @info "NaN in dTfdt" dTfdr d2Tfdr2 top_bound_term bot_bound_term Δξ[ir]
    end

    dTfdt[no_ice] .= 0 # Set to 0 elsewhere

    # lbound = [i for i in 2:dom.nr if (has_ice[i] && no_ice[i-1])]
    # rbound = [i for i in 1:dom.nr-1 if (has_ice[i] && no_ice[i+1])]
    # if length(lbound) == 0
    #     #nothing
    # elseif length(lbound) == 1
    #     # Treat ghost cell
    #     # dTfdt[lbound[1]-1] = ...
    # else
    #     @warn "Multiple chunks of ice may not be handled correctly." lbound rbound
    # end
    # if length(rbound) == 0
    #     #nothing
    # elseif length(rbound) == 1
    #     # Treat ghost cell
    #     # dTfdt[rbound[1]+1] = ...
    # else
    #     @warn "Multiple chunks of ice may not be handled correctly." lbound rbound
    # end

end


# -------------------------- Heat transfer only functions


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
    p_last = integ_pars[3]
    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)
    ϕ, Tf, Tw = ϕ_T_from_u(u, dom)

    T = solve_T(u, Tf, dom, params)

    vf, dϕdx_all = compute_frontvel_heat(u, T, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    dTf .= 0
    dTw .= 0

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



function dudt_heatmass_dae!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    p_last = integ_pars[3]
    # Tf_last = integ_pars[4]
    controls = integ_pars[5]

    input_measurements!(params, t, controls)

    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)

    # dTf .= 0 # Just in case some weird stuff got left there

    ϕ, Tf, Tw = ϕ_T_from_u_view(u, dom)
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    if minimum(ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[:l])/length(params[:l])
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[:Mw]/params[:R]/minT)
        flux = (calc_psub(minT) - params[:p_ch])/dom.zmax*b_ave
        dϕ .= flux/params[:ρf]
        # dϕ .= 0
        # dTw .= 0
        @info "no ice" extrema(dϕ)
        return nothing
    end

    # p = solve_p(u, T, dom, params, p_last)

    # Tf = pseudosteady_Tf_T_p(u, dom, params, Tf_last, p_last)
    # Tf = pseudosteady_Tf(u, dom, params, Tf_last)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)

    # Store pseudosteady in u, for use in other functions
    # Tf .= Tfs
    # dTf = (Tf - Tfg)/60

    integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    integ_pars[4] .= Tf
    vf, dϕdx_all = compute_frontvel_mass(u, Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # Compute time derivatives for 
    # dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    # dTfdt_radial!(dTf, u, T, p, dϕdx_all, dom, params)
    Qgl = compute_Qgl(u, T, dom, params)
    dTw .= (Q_gl_RF - Qgl) / m_cp_gl

    # dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all
    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        # dϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        dϕ[ind] = -rcomp - zcomp
    end

    dTfdt_radial!(dTf, u, Tf, T, p, dϕdx_all, dom, params)

    dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
    @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) extrema(T) Tw[1] params[:Tsh]
    if minimum(dϕ) < 0
        @info "negative dϕ" spy(ϕ .< 0) spy(dϕ .< 0) Tf[end]-Tf[1]
        # pl1 = heat(vr, dom)
        # pl2 = heat(vz, dom)
        # display(plot(pl1, pl2))
    end
    # if maximum(dϕ) > 1
    #     @info extrema(vz) extrema(vr) extrema(T) extrema(p) extrema(Tf)
    #     display(heat(T, dom))
    # end
    return nothing
end