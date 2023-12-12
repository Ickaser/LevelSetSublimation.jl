export dudt_heatmass, dudt_heatonly
export dudt_heatmass!, dudt_heatmass_dae!, dudt_heatmass_implicit!

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
    verbose = integ_pars[6]

    input_measurements!(params, t, controls)

    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)

    dTf .= 0 # Just in case some weird stuff got left there

    ϕ, Tw = ϕ_T_from_u_view(u, dom)[[true, false, true]]
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF, vial_thick = params

    if minimum(ϕ) > 0 # No ice left
        l_ave = 1/sum(1.0 ./params[:l])/length(params[:l])
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[:Mw]/params[:R]/minT)
        flux = (calc_psub(minT) - params[:p_ch])/dom.zmax*b_ave
        dϕ .= flux/params[:ρf]
        # dϕ .= 0
        # dTw .= 0
        verbose && @info "no ice, t=$t" extrema(dϕ)
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
    # A_vsh = π*((dom.rmax+vial_thick)^2 - dom.rmax^2)
    # Q_vsh = params[:Kv]*(params[:Tsh]-Tw[1]) * A_vsh
    # dTw .= (Q_gl_RF + Q_vsh - Qgl) / m_cp_gl
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
    if verbose && eltype(u) <: Float64
        dryfrac = 1 - compute_icevol_H(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) extrema(T) Tw[1] params[:Tsh]
        if minimum(dϕ) < 0
            # @info "negative dϕ" spy(ϕ .< 0) spy(dϕ .< 0) Tf[end]-Tf[1]
            # pl1 = heat(vr, dom)
            # pl2 = heat(vz, dom)
            # display(plot(pl1, pl2))
        end
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

    # isnan(qflux + ΔH*mflux) && @info "subheat" ir,iz qflux mflux dTdr dTdz dpdr dpdz

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
            # surf_area += compute_local_δ(cell, ϕ, dom)*locvol
        else # Cell is in ice, so need to extrapolate to get qloc
            # Identify possible neighbors across interface, take sum
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
            elseif length(nbs) == 1 # No vertical neighbor- just horizontal
                qloc = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nbs[1])..., dϕdx_all, dom, params)[1]
            elseif length(nbs) == 2 # One vertical neighbor, one horizontal
                qloc = mapreduce(+, nbs) do nb 
                    ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
                    Tuple(nb)[2] == 0 ?  ql*dϕdr : ql*dϕdz
                end
            elseif length(nbs) == 3 # Two vertical neighbor
                qloc = mapreduce(+, nbs) do nb 
                    ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
                    Tuple(nb)[2] == 0 ? ql*dϕdr : ql*dϕdz/2
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
                    sumfluxes += toparea * ΔH*(p_ch - calc_psub(Tf_loc))/Rp0 # Sublimation
                end
                if bot_contact[ir]
                    θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir+1,1])
                    ro = dom.rgrid[ir] + θr*dom.dr
                    botarea = π*(ro^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])

                    K_eff = 1/(1/params[:Kv] + Δξ[ir]/kf)
                    sumfluxes += botarea * K_eff*(params[:Tsh] - Tf_loc) # Shelf heat 
                end

                sumfluxes += Q_ic*vol # No A_l*kf*dTfdr becuase dTfdr=0
                # Write the energy balance differently for this case
                dTfdt[ir] = sumfluxes/vol/ρf/Cpf
                continue
            end
        elseif ir == dom.nr
            dTfdr = params[:Kw]/kf*(Tw - Tf[ir])
            if has_ice[dom.nr-1]
                d2Tfdr2 = (-2Tf[dom.nr] + 2Tf[dom.nr-1] + 2*dom.dr*dTfdr)*dom.dr2 # Robin ghost cell
            else
                # d2Tfdr2 = 0 
                # @warn "Possible mistreatment: set d2Tf/dr2 to 0 at right boundary, in an unlikely case"
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
                    sumfluxes += toparea * ΔH*(p_ch - calc_psub(Tf_loc))/Rp0 # Sublimation
                end
                if bot_contact[ir]
                    θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir-1,1])
                    ro = dom.rgrid[ir] 
                    ri = dom.rgrid[ir] - θr*dom.dr
                    botarea = π*(ro^2-ri^2)
                    Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                    K_eff = 1/(1/params[:Kv] + Δξ[ir]/kf)
                    sumfluxes += botarea * K_eff*(params[:Tsh] - Tf_loc) # Shelf heat 
                end

                A_l = dom.rgrid[ir] * Δξ[ir] *2π
                sumfluxes += -A_l*kf*dTfdr + Q_ic*vol
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
            # Q_surf_pp = 0.0
            # vol = 0.0
            # surf_area = 0.0
            # for cell in integ_cells
            #     locvol = dom.rgrid[Tuple(cell)[1]] * dom.dr * dom.dz * 2π
            #     if ϕ[cell] > 0
            #         qloc, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(cell)..., dϕdx_all, dom, params)
            #         # surf_area += compute_local_δ(cell, ϕ, dom)*locvol
            #     else # Cell is in ice, so need to extrapolate to get qloc
            #         # Identify possible neighbors across interface, take sum
            #         possible_nbs = [cell] .+ [CI(-1, 0), CI(0, 1), CI(0, -1)]
            #         nbs = [nb for nb in possible_nbs if (checkbounds(Bool, ϕ, nb) && ϕ[nb]>0)]
            #         if length(nbs) == 0
            #             @warn "um... no neighbors somehow?"
            #         elseif length(nbs) == 1 # No vertical neighbor- just horizontal
            #             qloc = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nbs[1])..., dϕdx_all, dom, params)[1]
            #         elseif length(nbs) == 2 # One vertical neighbor, one horizontal
            #             qloc = mapreduce(+, nbs) do nb 
            #                 ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
            #                 Tuple(nb)[2] == 0 ?  ql*dϕdr : ql*dϕdz
            #             end
            #         elseif length(nbs) == 3 # Two vertical neighbor
            #             qloc = mapreduce(+, nbs) do nb 
            #                 ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
            #                 Tuple(nb)[2] == 0 ? ql*dϕdr : ql*dϕdz/2
            #             end
            #         end
            #     end
            #     Q_surf_pp += qloc*compute_local_δ(cell, ϕ, dom)*locvol

            #     surf_area += compute_local_δ(cell, ϕ, dom)*locvol
            #     vol += compute_local_H(cell, ϕ, dom)*locvol
            # end
            sumfluxes = Q_surf_pp
            if top_contact[ir]
                θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir-1,dom.nz])
                ro = dom.rgrid[ir] + 0.5*dom.dr
                ri = dom.rgrid[ir] - θr*dom.dr
                toparea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                sumfluxes += toparea * ΔH*(p_ch - calc_psub(Tf_loc))/Rp0 # Sublimation
            end
            if bot_contact[ir]
                θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir-1,1])
                ro = dom.rgrid[ir] + 0.5*dom.dr
                ri = dom.rgrid[ir] - θr*dom.dr
                botarea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])

                K_eff = 1/(1/params[:Kv] + Δξ[ir]/kf)
                sumfluxes += botarea * K_eff*(params[:Tsh] - Tf_loc) # Shelf heat 
            end

            A_l = (dom.rgrid[ir] + 0.5dom.dr) * (Δξ[ir+1] + Δξ[ir])/2 *2π
            dTfdr = (Tf[ir+1] - Tf[ir])*dom.dr1 # We want derivative between grid points, so this is actually 2nd-order accurate
            sumfluxes +=  A_l*kf*dTfdr + Q_ic*vol
            # Write the energy balance differently for this case
            dTfdt[ir] = sumfluxes/vol/ρf/Cpf
            continue
        elseif no_ice[ir+1] # On right side: pulled away from wall
            topz = min(findlast(ϕ[ir,:] .< 0) + 1, dom.nz)
            botz = max(findfirst(ϕ[ir,:] .< 0)- 1, 1)
            integ_cells = [CI(iir, iz) for iz in botz:topz, iir in ir:ir+1]
            Q_surf_pp, surf_area, vol = Q_surf_integration(:right, integ_cells, ϕ, u, Tf, T, p, dϕdx_all, dom, params)
            # Q_surf_pp = 0.0
            # vol = 0.0
            # surf_area = 0.0
            # for cell in integ_cells
            #     locvol = dom.rgrid[Tuple(cell)[1]] * dom.dr * dom.dz * 2π
            #     if ϕ[cell] > 0
            #         qloc, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(cell)..., dϕdx_all, dom, params)
            #         # surf_area += compute_local_δ(cell, ϕ, dom)*locvol
            #     else # Cell is in ice, so need to extrapolate to get qloc
            #         # Identify possible neighbors across interface, take sum
            #         possible_nbs = [cell] .+ [CI(1, 0), CI(0, 1), CI(0, -1)]
            #         nbs = [nb for nb in possible_nbs if (checkbounds(Bool, ϕ, nb) && ϕ[nb]>0)]
            #         if length(nbs) == 0
            #             @warn "um... no neighbors somehow?"
            #         elseif length(nbs) == 1 # No vertical neighbor- just horizontal
            #             qloc = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nbs[1])..., dϕdx_all, dom, params)[1]
            #         elseif length(nbs) == 2 # One vertical neighbor, one horizontal
            #             qloc = mapreduce(+, nbs) do nb 
            #                 ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
            #                 Tuple(nb)[2] == 0 ?  ql*dϕdr : ql*dϕdz
            #             end
            #         elseif length(nbs) == 3 # Two vertical neighbor
            #             qloc = mapreduce(+, nbs) do nb 
            #                 ql, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, Tuple(nb)..., dϕdx_all, dom, params)
            #                 Tuple(nb)[2] == 0 ? ql*dϕdr : ql*dϕdz/2
            #             end
            #         else
            #             @warn "4 neighbors somehow in local surface integration"
            #         end
            #     end
            #     Q_surf_pp += qloc*compute_local_δ(cell, ϕ, dom)*locvol

            #     surf_area += compute_local_δ(cell, ϕ, dom)*locvol
            #     vol += compute_local_H(cell, ϕ, dom)*locvol
            # end
            sumfluxes = Q_surf_pp
            if top_contact[ir]
                θr = ϕ[ir,dom.nz] / (ϕ[ir,dom.nz] - ϕ[ir+1,dom.nz])
                ro = dom.rgrid[ir] + θr*dom.dr
                ri = dom.rgrid[ir] - 0.5*dom.dr
                toparea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                sumfluxes += toparea * ΔH*(p_ch - calc_psub(Tf_loc))/Rp0 # Sublimation
            end
            if bot_contact[ir]
                θr = ϕ[ir,1] / (ϕ[ir,1] - ϕ[ir+1,1])
                ro = dom.rgrid[ir] + θr*dom.dr
                ri = dom.rgrid[ir] - 0.5*dom.dr
                botarea = π*(ro^2-ri^2)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])

                K_eff = 1/(1/params[:Kv] + Δξ[ir]/kf)
                sumfluxes += botarea * K_eff*(params[:Tsh] - Tf_loc) # Shelf heat 
            end

            A_l = (dom.rgrid[ir] - 0.5dom.dr) * (Δξ[ir-1] + Δξ[ir])/2 *2π
            dTfdr = (Tf[ir] - Tf[ir-1])*dom.dr1 # We want derivative between grid points, so this is actually 2nd-order accurate
            sumfluxes += -A_l*kf*dTfdr + Q_ic*vol
            # Write the energy balance differently for this case
            dTfdt[ir] = sumfluxes/vol/ρf/Cpf
            continue
        else
            dTfdr = (Tf[ir+1] - Tf[ir-1])*0.5*dom.dr1
            d2Tfdr2 = (Tf[ir+1] - 2Tf[ir] + Tf[ir-1])*dom.dr2
        end

        # Treat top and bottom
        if top_contact[ir]
            # # adiabatic boundary 
            # top_bound_term = 0
            # direct sublimation boundary
            # top_bound_term = ΔH*(p_ch - calc_psub(Tf[ir]))/Rp0 * compute_local_H(CI(ir, dom.nz), ϕ, dom)*2
            top_bound_term = ΔH*(p_ch - calc_psub(Tf[ir]))/Rp0
            if typeof(top_bound_term) <: AbstractFloat
            end
        else
            iz = findlast(ϕ[ir,:] .<=0) + 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            top_bound_term = q - kf*dϕdr/dϕdz*dTfdr
            isnan(top_bound_term) && @info "top bound" q dϕdr dϕdz dTfdr
            # top_bound_term = q*dϕdz
        end

        if bot_contact[ir]
            # Shelf boundary
            # bot_bound_term = params[:Kv]*(Tf[ir] - params[:Tsh]) * compute_local_H(CI(ir, 1), ϕ, dom)*2
            K_eff = 1/(1/params[:Kv] + Δξ[ir]/kf)
            bot_bound_term = K_eff*(Tf[ir] - params[:Tsh]) 
            isnan(bot_bound_term) && @info "NaN bottom" Tf[ir] T[ir,begin] params[:Tsh]
        else
            # Stefan boundary
            iz = findfirst(ϕ[ir,:] .<=0) - 1
            q, dϕdr, dϕdz = local_sub_heating_dϕdx(u, Tf, T, p, ir, iz, dϕdx_all, dom, params)
            bot_bound_term = q - kf*dϕdr/dϕdz*dTfdr
            # bot_bound_term = q*dϕdz
            isnan(bot_bound_term) && @info "NaN bottom" q dϕdz dϕdr dTfdr
        end



        r1 = (ir == 1 ? 0 : 1/dom.rgrid[ir])
        dTfdt[ir] = (kf*(r1*dTfdr + d2Tfdr2) + Q_ic + 
            (top_bound_term - bot_bound_term)/Δξ[ir])/ρf/Cpf

        # isnan(dTfdt[ir]) && @info "NaN in dTfdt" dTfdr d2Tfdr2 top_bound_term bot_bound_term Δξ[ir]
    end

    dTfdt[no_ice] .= 0 # Set to 0 elsewhere

    # if any(isnan.(dTfdt))
    #     @info "NaN in dTfdt"
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



function dudt_heatmass_dae!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    # p_last = integ_pars[3]
    # Tf_last = integ_pars[4]
    controls = integ_pars[5]
    verbose = integ_pars[6]

    input_measurements!(params, t, controls)

    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)

    ϕ, Tf, Tw = ϕ_T_from_u_view(u, dom)
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    if any(Tf .< 0)
        @info "Negative temperatures passed"
        du .= NaN
        return
    end

    if minimum(ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[:l])/length(params[:l])
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[:Mw]/params[:R]/minT)
        flux = (calc_psub(minT) - params[:p_ch])/dom.zmax*b_ave
        dϕ .= flux/params[:ρf]
        # dϕ .= 0
        # dTw .= 0
        verbose && @info "no ice" extrema(dϕ)
        return nothing
    end


    # Tf = pseudosteady_Tf(u, dom, params, Tf_last)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)



    # integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    # integ_pars[4] .= Tf
    vf, dϕdx_all = compute_frontvel_mass(u, Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # Compute time derivatives for 
    # dTf .= Qice / ρf / Cpf / max(compute_icevol(ϕ, dom), 1e-6) # Prevent explosion during last time step by not letting volume go to 0
    # dTfdt_radial!(dTf, u, T, p, dϕdx_all, dom, params)
    dTfdt_radial!(dTf, u, Tf, T, p, dϕdx_all, dom, params)

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


    if verbose &&  eltype(u) <: Float64
        dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) extrema(T) Tw[1] params[:Tsh]
    end
    return nothing
end

function dudt_heatmass_implicit!(du, u, integ_pars, t)
    dom = integ_pars[1]
    params = integ_pars[2]
    # p_last = integ_pars[3]
    # Tf_last = integ_pars[4]
    controls = integ_pars[5]
    verbose = integ_pars[6]

    input_measurements!(params, t, controls)

    dϕ, dTf, dTw = ϕ_T_from_u_view(du, dom)

    ϕ, Tf, Tw = ϕ_T_from_u_view(u, dom)
    @unpack ρf, Cpf, m_cp_gl, Q_gl_RF = params

    if any(Tf .< 0)
        @info "Negative temperatures passed"
        du .= NaN
        return
    end

    if minimum(ϕ) > 0 # No ice left
        l_ave = 1/sum(1/params[:l])/length(params[:l])
        minT = minimum(solve_T(u, fill(NaN, dom.nr), dom, params))
        b_ave = l_ave * sqrt(params[:Mw]/params[:R]/minT)
        flux = (calc_psub(minT) - params[:p_ch])/dom.zmax*b_ave
        dϕ .= flux/params[:ρf]
        # dϕ .= 0
        # dTw .= 0
        verbose && @info "no ice" extrema(dϕ)
        return nothing
    end


    # Tf = pseudosteady_Tf(u, dom, params, Tf_last)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)



    # integ_pars[3] .= p # Cache current state of p as a guess for next timestep
    # integ_pars[4] .= Tf
    vf, dϕdx_all = compute_frontvel_mass(u, Tf, T, p, dom, params)
    extrap_v_fastmarch!(vf, u, dom)
    vr = @view vf[:, :, 1]
    vz = @view vf[:, :, 2]

    # Compute time derivatives for Tf
    dTfdt_radial!(dTf, u, Tf, T, p, dϕdx_all, dom, params)

    Qgl = compute_Qgl(u, T, dom, params)
    dTw .= (Q_gl_RF - Qgl) / m_cp_gl

    for ind in CartesianIndices(ϕ)
        ir, iz = Tuple(ind)

        dϕdr, dϕdz = choose_dϕdx_boundary(ir, iz, vr[ind] > 0, vz[ind] >0, dϕdx_all, dom)

        rcomp = dϕdr * vr[ind]
        zcomp = dϕdz * vz[ind]
        dϕ[ind] = -rcomp - zcomp
    end


    if verbose &&  eltype(u) <: Float64
        dryfrac = 1 - compute_icevol(ϕ, dom) / ( π* dom.rmax^2 *dom.zmax)
        @info "prog: t=$t, dryfrac=$dryfrac" extrema(dϕ) extrema(Tf) extrema(T) Tw[1] params[:Tsh]
    end
    return nothing
end