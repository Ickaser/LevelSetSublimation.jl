export dudt_heatmass, dudt_heatonly
export dudt_heatmass!, dudt_heatmass_dae!, dudt_heatmass_implicit!

"""
    dudt_heatmass!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

Splitting `u` and `du` into `ϕ`, `Tf`, and `Tvw` is handled by `ϕ_T_from_u` and `ϕ_T_from_u_view`.

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function dudt_heatmass!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params_vary = integ_pars[2]
    verbose = integ_pars[3]
    saved_Tf = integ_pars[4]

    params = params_vary(t)

    du.Tf .= 0 # Just in case some weird stuff got left there

    @unpack ρf, Cpf, ρ_vw, cp_vw = params[1]
    @unpack A_v, m_v = params[2]
    QRFvw = calc_QpppRFvw(params) * m_v/ρ_vw

    if minimum(u.ϕ) > 0 # No ice left
        l_ave = 1/sum(1.0 ./params[2].l)/length(params[2].l)
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[1].Mw/params[1].R/minT)
        flux = (calc_psub(minT) - params[3].pch)/dom.zmax*b_ave
        du.ϕ .= flux/params[1].ρf
        du.Tvw = flux*params[1].ΔH*dom.rmax^2*π
        verbose && @info "no ice, t=$t" extrema(du.ϕ)
        return nothing
    end

    # Tf_g = Tf_last
    # Tf_g = u.Tf
    Tf_g = Tf_guess(u.Tf, t, saved_Tf)
    Tf = pseudosteady_Tf(u, dom, params, Tf_g)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)

    # Save the computed Tf here, as well as in the callback
    # integ_pars[4] .= Tf

    vf, dϕdx_all = compute_frontvel_mass(u, Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    for ind in CartesianIndices(u.ϕ)
        ir, iz = Tuple(ind)
        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)
        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        du.ϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
    end

    # Compute time derivative for Tvw
    Q_vwf = compute_Qvwf(u, T, dom, params)
    A_p = π*dom.rmax^2
    Q_shvw = (A_v-A_p) * params[3].Kshf * (params[3].Tsh - u.Tvw)
    du.Tvw = (QRFvw - Q_vwf + Q_shvw) / m_v / cp_vw

    if verbose && eltype(u) <: Float64
        dryfrac = 1 - compute_icevol_H(u.ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(du.ϕ) extrema(Tf) extrema(T) u.Tvw params[3].Tsh
    end
    return nothing
end

"""
doc
"""
function dudt_heatmass_dae!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params_vary = integ_pars[2]
    verbose = integ_pars[3]

    params = params_vary(t)

    @unpack ρf, Cpf, ρ_vw, cp_vw = params[1]
    @unpack A_v, m_v = params[2]
    QRFvw = calc_QpppRFvw(params) * m_v/ρ_vw

    if any(u.Tf .< 0)
        verbose && @info "Negative temperatures passed"
        du .= NaN
        return
    end

    if minimum(u.ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[2].l)/length(params[2].l)
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[1].Mw/params[1].R/minT)
        flux = (calc_psub(minT) - params[3].pch)/dom.zmax*b_ave
        du.ϕ .= flux/params[1].ρf
        # du.ϕ .= 0
        # du.Tvw .= 0
        verbose && @info "no ice" extrema(du.ϕ)
        return nothing
    end


    T = solve_T(u, u.Tf, dom, params)
    p = solve_p(u, u.Tf, T, dom, params)

    vf, dϕdx_all = compute_frontvel_mass(u, u.Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    for ind in CartesianIndices(u.ϕ)
        ir, iz = Tuple(ind)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        du.ϕ[ind] = max(0.0, -rcomp - zcomp) # Prevent solidification
        # du.ϕ[ind] = -rcomp - zcomp
    end

    # du.Tf .= 0 # zero out in case there was junk there before
    dTfdt_radial!(du.Tf, u, u.Tf, T, p, dϕdx_all, dom, params)

    Q_vwf = compute_Qvwf(u, T, dom, params)
    A_p = π*dom.rmax^2
    Q_shvw = (A_v-A_p) * params[3].Kshf * (params[3].Tsh - u.Tvw)
    du.Tvw = (QRFvw - Q_vwf + Q_shvw) / m_v / cp_vw



    if verbose &&  eltype(u) <: Float64
        dryfrac = 1 - compute_icevol_H(u.ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        Δξ = compute_iceht_bottopcont(u.ϕ, dom)[1]
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(du.ϕ) extrema(u.Tf) extrema(T) u.Tvw params[3].Tsh extrema(Δξ) extrema(du.Tf)
    end
    return nothing
end

"""
doc
"""
function dudt_heatmass_implicit!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params_vary = integ_pars[2]
    verbose = integ_pars[3]

    params = params_vary(t)

    @unpack ρf, Cpf, ρ_vw, cp_vw = params[1]
    @unpack A_v, m_v = params[2]
    QRFvw = calc_QpppRFvw(params) * m_v/ρ_vw

    if any(u.Tf .< 0)
        verbose && @info "Negative temperatures passed"
        du .= NaN
        return
    end

    if minimum(u.ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[2].l)/length(params[2].l)
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[1].Mw/params[1].R/minT)
        flux = (calc_psub(minT) - params[3].pch)/dom.zmax*b_ave
        du.ϕ .= flux/params[1].ρf
        # dϕ .= 0
        # du.Tvw .= 0
        verbose && @info "no ice" extrema(du.ϕ)
        return nothing
    end


    T = solve_T(u, u.Tf, dom, params)
    p = solve_p(u, u.Tf, T, dom, params)



    vf, dϕdx_all = compute_frontvel_mass(u, u.Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # Compute time derivatives for Tf
    dTfdt_radial!(du.Tf, u, u.Tf, T, p, dϕdx_all, dom, params)

    Q_vwf = compute_Qvwf(u, T, dom, params)
    A_p = π*dom.rmax^2
    Q_shvw = (A_v-A_p) * params[3].Kshf * (params[3].Tsh - u.Tvw)
    du.Tvw .= (QRFvw - Q_vwf + Q_shvw) / m_v / cp_vw

    for ind in CartesianIndices(u.ϕ)
        ir, iz = Tuple(ind)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        du.ϕ[ind] = max(0.0, -rcomp - zcomp)
        # du.ϕ[ind] = -rcomp - zcomp
    end


    if verbose &&  eltype(u) <: Float64
        dryfrac = 1 - compute_icevol_H(u.ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(du.ϕ) extrema(u.Tf) extrema(T) u.Tvw params[3].Tsh
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
    # return dϕdr, dϕdz
    normed = hypot(dϕdr, dϕdz)
    return dϕdr/normed, dϕdz/normed
end


"""
    dudt_heatmass(u, dom::Domain, params, t)
    dudt_heatmass(u, config::Dict, t=0.0)
    
Compute the time derivative of `u` with given (nondimensional) parameters.

Wraps a call on `dudt_heatmass!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatmass(u, dom::Domain, params_n, t)
    integ_pars = (dom, params_n, true, similar(u, dom.nr))
    du = similar(u)
    dudt_heatmass!(du, u, integ_pars, t)
    return du
end
function dudt_heatmass(u, config::Dict, t=0.0)
    # Set up simulation domain & parameters
    dom = Domain(config)
    params_vary = params_nondim_setup(config[:cparams], config[:controls])
    params = params_vary(t)
    dudt_heatmass(u, dom, params)
end

"""
    dudt_heatmass_params(u, config)
    
Compute the time derivative of `u` with given parameters, and also return `dom` and `params` associated with the given `config`.

Wraps a call on `dudt_heatmass`, for convenience in debugging and elsewhere that efficiency is less important
"""
function dudt_heatmass_params(u, config, t)
    # Set up simulation domain & parameters
    dom = Domain(config)
    paramsn = params_nondim_setup(config[:cparams], config[:controls])

    dudt_heatmass(u, dom, paramsn, t), dom, paramsn
end

# ---------------------------
# For pseudosteady radial temperature

function local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
    @unpack ΔH = params[1]
    @unpack kd = params[2]
    b = eval_b_loc(T, p, ir, iz, params)

    dTdr, dTdz = compute_Tderiv(u, Tf, T, ir, iz, dom, params)
    dpdr, dpdz = compute_pderiv(u, Tf, T, p, ir, iz, dom, params)
    dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, dpdr<0, dpdz<0, dϕdx_all, dom)

    qflux = kd*(dϕdr*dTdr + dϕdz*dTdz)
    # mflux = b*(dϕdr*dpdr + dϕdz*dpdz)
    mflux = min(0.0, b*(dϕdr*dpdr + dϕdz*dpdz)) # Prevent desublimation in energy balance
    return qflux + ΔH*mflux, dϕdr, dϕdz
end


function Q_surf_integration(side, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
    Q_surf_pp = 0.0
    vol = 0.0
    surf_area = 0.0
    for cell in integ_cells
        ir = Tuple(cell)[1]
        # locvol = dom.rgrid[ir] * dom.dr * dom.dz * 2π
        # if locvol == 0 && ir == 1
        #     locvol = (0.5dom.dr)^2 *2π * dom.dz
        # end
        ri = max(0, dom.rgrid[ir] -0.5dom.dr)
        ro = min(dom.rmax, dom.rgrid[ir] +0.5dom.dr)
        locvol = (ro^2-ri^2) * dom.dz * π
        if ϕ[cell] > 0
            qloc, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(cell)..., dϕdx_all, dom, params)
        else # Cell is in ice, so need to extrapolate to get qloc
            # Identify possible neighbors across interface, take sum
            # At most one side neighbor, in normal circumstances
            if side == :left
                possible_nbs = [cell] .+ [CI(0, 1), CI(0, -1), CI(-1, 0)]
            elseif side == :right
                possible_nbs = [cell] .+ [CI(1, 0), CI(0, 1), CI(0, -1)]
            else
                @warn "internal argument passed weird: :left or :right" side
                possible_nbs = [cell] .+ [CI(1, 0), CI(0, 1), CI(0, -1), CI(-1, 0)]
            end
            nbs = [nb for nb in possible_nbs if (checkbounds(Bool, ϕ, nb) && ϕ[nb]>0)]
            if length(nbs) == 0
                @warn "um... no neighbors somehow?"
            elseif length(nbs) == 1 # One neighbor: just use its value
                qloc = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nbs[1])..., dϕdx_all, dom, params)[1]
            elseif length(nbs) == 2 # At most one horizontal neighbor, so probably one vertical and one horizontal: take a dot product
                qloc = mapreduce(+, nbs) do nb 
                    ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
                    Tuple(nb)[2] == Tuple(cell)[2] ?  ql*dϕdr : ql*dϕdz
                end
            elseif length(nbs) == 3 # Three neighbors guarantees one horizontal and two vertical
                qloc = mapreduce(+, nbs) do nb 
                    ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
                    Tuple(nb)[2] == Tuple(cell[2]) ? ql*dϕdr : ql*dϕdz/2
                end
            else
                @warn "4 neighbors somehow in local surface integration"
            end
        end
        Q_surf_pp += qloc*compute_local_δ(cell, ϕ, dom)*locvol

        surf_area += compute_local_δ(cell, ϕ, dom)*locvol
        vol += compute_local_H(cell, ϕ, dom)*locvol
    end
    return Q_surf_pp, surf_area, vol
end

function dTfdt_radial(u, Tf, T, p, dϕdx_all, dom::Domain, params)
    dTfdt = similar(u, dom.nr)
    dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom, params)
    return dTfdt
end

function dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom::Domain, params)
    ϕ = u.ϕ
    Tvw = u.Tvw
    @unpack ρf, Cpf, kf, ΔH = params[1]
    @unpack Rp0, Kvwf = params[2]
    @unpack Tsh, pch, Kshf = params[3]
    QRFf = calc_QpppRFf.(Tf, [params])

    Δξ, bot_contact, top_contact = compute_iceht_bottopcont(ϕ, dom)

    has_ice = (Δξ .> 0)
    no_ice = (Δξ .== 0)

    for ir in axes(dTfdt, 1)
        if no_ice[ir]
            continue
        end

        # Treat radial derivatives

        if ir == 1
            dTfdr = 0
            if has_ice[2]
                d2Tfdr2 = (-2Tf[1] + 2Tf[2])*dom.dr2 # Adiabatic ghost cell
            else
                # d2Tfdr2 = 0 
                # @warn "Possible mistreatment: set d2Tf/dr2 to 0 at left boundary, in an unlikely case"
                topz = min(findlast(ϕ[ir,:] .< 0) + 1, dom.nz)
                botz = max(findfirst(ϕ[ir,:] .< 0)- 1, 1)
                integ_cells = [CI(iir, iz) for iz in botz:topz, iir in ir:ir+1]
                Q_surf_pp, surf_area, vol = Q_surf_integration(:right, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
                sumfluxes = Q_surf_pp
                if top_contact[ir]
                    θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir+1,dom.nz])
                    ro = dom.rgrid[ir] + θr*dom.dr
                    toparea = π*(ro^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                    sumfluxes += toparea * ΔH*(pch - calc_psub(Tf_loc))/Rp0 # Sublimation
                end
                if bot_contact[ir]
                    θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir+1,1])
                    ro = dom.rgrid[ir] + θr*dom.dr
                    botarea = π*(ro^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])

                    # K_eff = 1/(1/Kshf + Δξ[ir]/kf)
                    K_eff = Kshf
                    sumfluxes += botarea * K_eff*(Tsh - Tf_loc) # Shelf heat 
                end

                sumfluxes += QRFf[ir]*vol # No A_l*kf*dTfdr becuase dTfdr=0
                # Write the energy balance differently for this case
                dTfdt[ir] = sumfluxes/vol/ρf/Cpf
                continue
            end
        elseif ir == dom.nr
            dTfdr = Kvwf/kf*(Tvw - Tf[ir])
            # if Δξ[ir] < dom.dz/2 # At small ice height, need a special case
            if Δξ[ir]/dom.dz < θ_THRESH # At small ice height, need a special case
                iz = findlast(ϕ[ir,:] .<=0) + 1
                q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)

                qtop = q/dϕdz

                # We get terms that go to zero, and can therefore solve a linear equaiton for Tf, other than q_top
                # den = (Δξ[ir]/dom.dr/2*Kvwf + Kshf)
                # rhs = (Δξ[ir]/dom.dr/2*Kvwf*Tvw + Kshf*Tsh + qtop)/den
                den = Kshf
                rhs = (Kshf*Tsh + qtop)/den
                # dTfdt[ir] = rhs - Tf[ir]
                # timescale = Δξ[ir]^2/kf*ρf*Cpf
                timescale = dom.zmax^2/kf*ρf*Cpf
                dTfdt[ir] = (rhs - Tf[ir])/timescale
                # @info "Small ice at corner" iz ir Δξ[ir] timescale dTfdt[ir-3:ir] rhs Tf[ir]
                continue

            end
            # elseif has_ice[dom.nr-1]
            if has_ice[dom.nr-1]
                d2Tfdr2 = (-2Tf[dom.nr] + 2Tf[dom.nr-1] + 2*dom.dr*dTfdr)*dom.dr2 # Robin ghost cell
            else
                topz = min(findlast(ϕ[ir,:] .< 0) + 1, dom.nz)
                botz = max(findfirst(ϕ[ir,:] .< 0)- 1, 1)
                integ_cells = [CI(iir, iz) for iz in botz:topz, iir in ir-1:ir]
                Q_surf_pp, surf_area, vol = Q_surf_integration(:left, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
                sumfluxes = Q_surf_pp
                if top_contact[ir]
                    θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir-1,dom.nz])
                    ro = dom.rgrid[ir]
                    ri = dom.rgrid[ir] - θr*dom.dr
                    toparea = π*(ro^2-ri^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                    sumfluxes += toparea * ΔH*(pch - calc_psub(Tf_loc))/Rp0 # Sublimation
                end
                if bot_contact[ir]
                    θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir-1,1])
                    ro = dom.rgrid[ir] 
                    ri = dom.rgrid[ir] - θr*dom.dr
                    botarea = π*(ro^2-ri^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                    # K_eff = 1/(1/Kshf + Δξ[ir]/kf)
                    K_eff = Kshf
                    sumfluxes += botarea * K_eff*(Tsh - Tf_loc) # Shelf heat 
                end

                A_l = dom.rgrid[ir] * Δξ[ir] *2π
                sumfluxes += -A_l*kf*dTfdr + QRFf[ir]*vol
                # Write the energy balance differently for this case
                dTfdt[ir] = sumfluxes/vol/ρf/Cpf
                continue
            end
        elseif no_ice[ir-1] && no_ice[ir+1] # Ice on both sides
            # @warn "Ice surrounded by gap: ignore radial gradients, treat only vertical." has_ice[ir-1:ir+1] 
            dTfdr = 0
            d2Tfdr2 = 0
            # dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            # d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        elseif no_ice[ir-1] # On left side: away from center
            topz = min(findlast(ϕ[ir,:] .< 0) + 1, dom.nz)
            botz = max(findfirst(ϕ[ir,:] .< 0)- 1, 1)
            integ_cells = [CI(iir, iz) for iz in botz:topz, iir in ir-1:ir]
            Q_surf_pp, surf_area, vol = Q_surf_integration(:left, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
            sumfluxes = Q_surf_pp
            if top_contact[ir]
                θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir-1,dom.nz])
                ro = dom.rgrid[ir] + 0.5*dom.dr
                ri = dom.rgrid[ir] - θr*dom.dr
                toparea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                sumfluxes += toparea * ΔH*(pch - calc_psub(Tf_loc))/Rp0 # Sublimation
            end
            if bot_contact[ir]
                θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir-1,1])
                ro = dom.rgrid[ir] + 0.5*dom.dr
                ri = dom.rgrid[ir] - θr*dom.dr
                botarea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])

                # K_eff = 1/(1/Kshf + Δξ[ir]/kf)
                K_eff = Kshf
                sumfluxes += botarea * K_eff*(Tsh - Tf_loc) # Shelf heat 
            end

            A_l = (dom.rgrid[ir] + 0.5dom.dr) * (Δξ[ir+1] + Δξ[ir])/2 *2π
            dTfdr = (Tf[ir+1] - Tf[ir])*dom.dr1 # We want derivative between grid points, so this is actually 2nd-order accurate
            sumfluxes +=  A_l*kf*dTfdr + QRFf[ir]*vol
            # Write the energy balance differently for this case
            dTfdt[ir] = sumfluxes/vol/ρf/Cpf
            continue
        elseif no_ice[ir+1] # On right side: pulled away from wall
            topz = min(findlast(ϕ[ir,:] .< 0) + 1, dom.nz)
            botz = max(findfirst(ϕ[ir,:] .< 0)- 1, 1)
            integ_cells = [CI(iir, iz) for iz in botz:topz, iir in ir:ir+1]
            Q_surf_pp, surf_area, vol = Q_surf_integration(:right, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
            sumfluxes = Q_surf_pp
            if top_contact[ir]
                θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir+1,dom.nz])
                ro = dom.rgrid[ir] + θr*dom.dr
                ri = dom.rgrid[ir] - 0.5*dom.dr
                toparea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                sumfluxes += toparea * ΔH*(pch - calc_psub(Tf_loc))/Rp0 # Sublimation
            end
            if bot_contact[ir]
                θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir+1,1])
                ro = dom.rgrid[ir] + θr*dom.dr
                ri = dom.rgrid[ir] - 0.5*dom.dr
                botarea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])

                # K_eff = 1/(1/Kshf + Δξ[ir]/kf)
                K_eff = Kshf
                sumfluxes += botarea * K_eff*(Tsh - Tf_loc) # Shelf heat 
            end

            A_l = (dom.rgrid[ir] - 0.5dom.dr) * (Δξ[ir-1] + Δξ[ir])/2 *2π
            dTfdr = (Tf[ir] - Tf[ir-1])*dom.dr1 # We want derivative between grid points, so this is actually 2nd-order accurate
            sumfluxes += -A_l*kf*dTfdr + QRFf[ir]*vol
            # Write the energy balance differently for this case
            dTfdt[ir] = sumfluxes/vol/ρf/Cpf
            continue
        else
            dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        end

        # Treat top and bottom
        if top_contact[ir]
            # direct sublimation boundary
            # top_bound_term = ΔH*(pch - calc_psub(Tf[ir]))/Rp0 * compute_local_H(CI(ir, dom.nz), ϕ, dom)*2
            top_bound_term = ΔH*(pch - calc_psub(Tf[ir]))/Rp0
        else
            iz = findlast(ϕ[ir,:] .<=0) + 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            # top_bound_term = q - kf*dϕdr/dϕdz*dTfdr
            dξdr = -dϕdr/dϕdz
            top_bound_term = q/dϕdz + kf*dξdr*dTfdr
        end

        if bot_contact[ir]
            # Shelf boundary
            # bot_bound_term = Kshf*(Tf[ir] - Tsh) * compute_local_H(CI(ir, 1), ϕ, dom)*2
            # if Δξ[ir] < dom.dz
            #     K_eff = Kshf
            # else
            #     K_eff = 1/(1/Kshf + Δξ[ir]/kf)
            # end
            K_eff = Kshf
            bot_bound_term = K_eff*(Tf[ir] - Tsh) 
        else
            # Stefan boundary
            iz = findfirst(ϕ[ir,:] .<=0) - 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            dξdr = -dϕdr/dϕdz
            bot_bound_term = q/dϕdz + kf*dξdr*dTfdr
        end

        r1 = (ir == 1 ? 0 : 1/dom.rgrid[ir])
        dTfdt[ir] = (kf*(r1*dTfdr + d2Tfdr2) + QRFf[ir] + 
            (top_bound_term - bot_bound_term)/Δξ[ir])/ρf/Cpf
    end
    dTfdt[no_ice] .= 0 # Set to 0 elsewhere
end
