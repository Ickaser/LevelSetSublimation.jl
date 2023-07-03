export compute_Qice, compute_icesurf, compute_icevol
export compute_frontvel_heat, compute_frontvel_mass, plot_frontvel
export compute_frontvel_fixedspeed

"""
    function compute_Qice(ϕ, T, p, dom::Domain, params)

Compute the total heat input into frozen & dried domains. Also passes glass-ice heat as separate return.

See p. 112-113 from paperlike notes. At pseudosteady conditions, all heat getting added to the system goes to the frozen domain,
so we don't actually need to treat different areas of the surface.
In addition, the sublimation flux is simply evaluated on the top surface.

"""
function compute_Qice(u, T, p, dom::Domain, params)

    @unpack ΔH= params
    # ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)[1]

    Q_glshvol, Qgl = compute_Qice_noflow(u, T, dom, params)

    # Sublimation rate
    md = compute_topmassflux(u, T, p, dom, params)
    Qsub = - md * ΔH

    return Q_glshvol + Qsub, Qgl
end


"""
    compute_Qice_noflow(u, T, dom::Domain, params)

Compute the total heat input into frozen & dried domains, assuming mass flow is zero. Also passes glass-ice heat as separate return.
"""
function compute_Qice_noflow(u, T, dom::Domain, params)

    @unpack Kv, Kgl, Q_ic, Q_ck, Tsh= params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)

    # Heat flux from shelf, at bottom of vial
    rweights = ones(Float64, dom.nr) 
    rweights[begin] = rweights[end] = 0.5
    Qsh = 2π* Kv * dom.dr* sum(rweights .* dom.rgrid .* (Tsh .- T[:,1] ))

    # Heat flux from glass, at outer radius
    zweights = ones(Float64, dom.nz) 
    zweights[begin] = zweights[end] = 0.5
    Qgl = 2π*dom.rmax * Kgl * dom.dz * sum(zweights .* ( Tgl .- T[end,:]))

    # Volumetric heat in cake and ice
    icevol = compute_icevol(ϕ, dom)
    dryvol = π*dom.rmax^2*dom.zmax - icevol
    Qvol = icevol * Q_ic + dryvol * Q_ck
    
    # @info "Q" Tsh Tf Qsh Qgl Qvol icevol

    return Qsh + Qgl + Qvol, Qgl
end

"""
    compute_Qice_nodry(u, T, dom::Domain, params)

Compute the total heat input into frozen domain from volumetric, shelf, and glass. 
Contrast with `compute_Qice_noflow` and `compute_Qice`, which include heat to dried domain.
"""
function compute_Qice_nodry(u, T, dom::Domain, params)
    @unpack Kv, Kgl, Q_ic, Q_ck, Tsh = params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)

    # Heat flux from shelf, at bottom of vial
    rweights = compute_icesh_area_weights(ϕ, dom)
    Qsh = 2π* Kv * sum(rweights .* (Tsh .- T[:,1] ))
    

    # Heat flux from glass, at outer radius
    zweights = compute_icegl_area_weights(ϕ, dom)
    Qgl = 2π*dom.rmax * Kgl * sum(zweights .* ( Tgl .- T[end,:]))

    # Volumetric heat in cake and ice
    icevol = compute_icevol(ϕ, dom)
    Qvol = icevol * Q_ic 

    # @info "geom" icevol get_subf_r(ϕ, dom)
    return Qsh + Qgl + Qvol
end

function compute_Qgl(u, T, dom::Domain, params)
    @unpack Kgl = params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    # Heat flux from glass, at outer radius
    zweights = compute_icegl_area_weights(ϕ, dom)
    Qgl = 2π*dom.rmax * Kgl * sum(zweights .* ( Tgl .- T[end,:]))
end

"""
    compute_topmassflux(ϕ, T, p, dom::Domain, params)

Compute total mass flow through top of the cake (that is, mass flux integrated across top surface).
"""
function compute_topmassflux(u, T, p, dom::Domain, params)
    dpdz = [compute_pderiv(u, T, p, ir, dom.nz, dom, params)[2] for ir in 1:dom.nr]
    b = eval_b(T, p, params)[:,end] # all r, top of z
    rweights = zeros(Float64, dom.nr) 
    for ir in 1:dom.nr-1
        rmid = (dom.rgrid[ir] + dom.rgrid[ir+1])/2
        rweights[ir] += rmid^2/2
        rweights[ir+1] -= rmid^2/2
    end
    # md = - 2π * sum(rweights .* b .*dpdz )
    md = - 2π * sum(dpdz .*b .* dom.rgrid) * dom.dr
    # @info "massflux" dpdz b md 1/dom.dz
    return md
end

"""
    compute_Tderiv(u, T, ir::Int, iz::Int, dom::Domain, params)

Compute (`∂T/∂r`, `∂T/∂z`) at point `(ir, iz)`, given system state `u`, `T`.

Quadratic ghost cell extrapolation (into frozen domain), second order finite differences, for T.
"""
function compute_Tderiv(u, Tf, T, ir::Int, iz::Int, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, nr, nz = dom
    @unpack k, Kv, Kgl, Tsh = params

    ϕ, Tgl = ϕ_T_from_u(u, dom)[[true, false, true]]
    pT = T[ir, iz]
    ϕp = ϕ[ir, iz]
    r = dom.rgrid[ir]
    
    if ϕp > 2dr || ϕp > 2dz || ϕp < -2dr || ϕp < -2dz
        @debug "Computing heat flux for cell which may not be at front." ir iz ϕp
    end

    # θ_thresh = max(1/nr, 1/nz)
    θ_thresh = 0.05

    # Enforce BCs explicitly for boundary cells
    if ir == 1 # Symmetry
        dTr = 0
    elseif ir == nr # Robin: glass
        dTr = Kgl/k*(Tgl - pT)
    else 
        # Bulk
        eϕ = ϕ[ir+1, iz]
        wϕ = ϕ[ir-1, iz]

        eT = T[ir+1, iz]
        wT = T[ir-1, iz]
        if wϕ <= 0 && eϕ <= 0
            # Constant Tf:
            # θr1 = ϕp/(ϕp - wϕ)
            # θr2 = ϕp/(ϕp - eϕ)
            # dTr = 0.5*dr1*((Tf - pT)*(1/θr2 - 1/θr1)
            # Varying Tf: pretend in bulk
            dTr = 0.5*dr1*(Tf[ir+1]-Tf[ir-1])
        elseif wϕ <= 0 # West ghost cell
            θr = ϕp /(ϕp - wϕ)
            Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
            if θr > θ_thresh
                # dTr = (-Tf/(1+θr)/θr + pT*(1-θr)/θr + eT*(θr)/(θr+1)) * dr1 # Quadratic extrapolation
                dTr = (pT - Tf_loc)/θr * dr1 # LInear extrapolation
            else 
                dTr = (eT - Tf_loc)/(θr+1)*dr1 # Linear extrapolation from east
            end
        elseif eϕ <= 0 # East ghost cell
            θr = ϕp /(ϕp - eϕ)
            Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
            if θr > θ_thresh
                # dTr = ( Tf/(θr+1)/θr - pT*(1-θr)/θr - wT*θr/(θr+1)) * dr1 # Quadratic extrapolation
                dTr = (Tf_loc - pT)/θr * dr1 # LInear extrapolation
            else
                dTr = (Tf_loc - wT)/(θr+1)*dr1 # Linear extrapolation from west
            end
        else # No ghost cells
            dTr = (eT - wT) * 0.5*dr1 # Centered difference
        end
    end
            
    # For all z derivatives, use Tf[ir]
    Tf_loc = Tf[ir]
    # Enforce BCs explicitly for boundary cells
    if iz == 1 # Robin: shelf
        dTz = Kv/k*(Tsh-pT)
    elseif iz == nz # Adiabatic
        dTz = 0
    else 
        # Bulk
        nϕ = ϕ[ir, iz+1]
        sϕ = ϕ[ir, iz-1]

        nT = T[ir, iz+1]
        sT = T[ir, iz-1]
        # North and south ghost cell: weird kink
        if sϕ <= 0 && nϕ <= 0
            θz1 = ϕp/(ϕp - nϕ)
            θz2 = ϕp/(ϕp - sϕ)
            dTz = 0.5*dz1*(Tf_loc - pT)*(1/θz2 - 1/θz1)
        elseif sϕ <= 0 # South ghost cell
            θz = ϕp /(ϕp - sϕ)
            if θz > θ_thresh
                # dTz = (-Tf_loc/(θz+1)/θz + pT*(1-θz)/θz + nT*θz/(θz+1)) * dz1 # Quadratic extrapolation
                dTz = (pT - Tf_loc)/θz * dz1 # Linear extrapolation
            else
                dTz = (nT - Tf_loc)/(θz+1)*dz1
            end
        elseif nϕ <= 0 # North ghost cell
            θz = ϕp /(ϕp - nϕ)
            if θz > θ_thresh
                # dTz = ( Tf_loc/(θz+1)/θz - pT*(1-θz)/θz - sT*θz/(θz+1)) * dz1 # Quadratic extrapolation
                dTz = (Tf_loc - pT)/θz * dz1 # Linear extrapolation
            else
                dTz = (Tf_loc - sT )/(θz+1)*dz1
            end
        else # No ghost cells
            dTz = (nT - sT) * 0.5*dz1 # Centered difference
        end
    end

    return dTr, dTz
end

"""
    compute_pderiv(u, T, p, ir::Int, iz::Int, dom::Domain, params)

Compute  `∂p/∂r`, `∂p/∂z` at point `(ir, iz)`, given system state `u`, `T`, `p`.

Quadratic ghost cell extrapolation (into frozen domain), second order finite differences, for p.

The distinction in discretization between this and `compute_Tderiv` is essentially just the boundary treatments.
"""
function compute_pderiv(u, Tf, T, p, ir::Int, iz::Int, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, nr, nz = dom
    @unpack p_ch, Rp0 = params
    ϕ, Tgl = ϕ_T_from_u(u, dom)[[true, false, true]]
    pp = p[ir, iz]
    ϕp = ϕ[ir, iz]
    
    if ϕp > 2dr || ϕp > 2dz || ϕp < -2dr || ϕp < -2dz
        @debug "Computing mass flux for cell which may not be at front." ir iz ϕp
    end

    # θ_thresh = max(1/nr, 1/nz)
    θ_thresh = 0.05

    # Enforce BCs explicitly for boundary cells
    if ir == 1 
        dpr = 0
    elseif ir == nr
        dpr = 0
    else 
        # Bulk
        eϕ = ϕ[ir+1, iz]
        wϕ = ϕ[ir-1, iz]
        pp = p[ir, iz]
        pe = p[ir+1, iz]
        pw = p[ir-1, iz]
        dpr = (pe - pw) * 0.5dr1

        if wϕ <= 0 && eϕ <= 0 # West and east ghost cell
            # Constant Tf
            # θr1 = ϕp/(ϕp - wϕ)
            # θr2 = ϕp/(ϕp - eϕ)
            # dpr = 0.5*dr1*(psub_l - pp)*(1/θr2 - 1/θr1)
            # Varying Tf; I'm lazy and this is a rare case
            dpr = 0.5*dr1*(calc_psub(Tf[ir+1]) - calc_psub(Tf[ir-1]))
        elseif wϕ <= 0 # West ghost cell
            θr = ϕp /(ϕp - wϕ)
            Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
            psub_l = calc_psub(Tf_loc)
            if θr > θ_thresh
                dpr = (-psub_l/(1+θr)/θr + pp*(1-θr)/θr + pe*(θr)/(θr+1)) * dr1 # Quadratic extrapolation
                # dpr = (pp - psub_l)/θr*dr1 # Linear extrapolation
            else 
                dpr = (pe - psub_l)/(θr+1)*dr1 # Linear extrapolation, further out
            end
        elseif eϕ <= 0 # East ghost cell
            θr = ϕp /(ϕp - eϕ)
            Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
            psub_l = calc_psub(Tf_loc)
            if θr > θ_thresh
                dpr = ( psub_l/(θr+1)/θr - pp*(1-θr)/θr - pw*θr/(θr+1)) * dr1 # Quadratic extrapolation
                # dpr = (psub_l - pp)/θr*dr1 # Linear extrapolation
            else
                dpr = (psub_l - pw)/(θr+1)*dr1 # Linear extrapolation, further out
            end
        else # No ghost cells
            dpr = (pe - pw) * 0.5*dr1 # Centered difference
        end
        # if (ir, iz) == (43, 38)
        # @info "pderiv" ir iz wϕ ϕp eϕ  pw pp pe dpr
        # end

    end
            
    # All z derivatives use Tf[ir]
    psub_l = calc_psub(Tf[ir])
    if iz == 1 
        # Enforce BCs explicitly for Neumann boundary cells
        dpz = 0
    elseif iz == nz 
        # Robin boundary condition: employ explicitly
        bp = eval_b_loc(T, p, ir, iz, params)
        # bp*dpz = Δp/Rp0
        dpz = (p_ch - p[ir,iz])/bp/Rp0
        # @info "here" ir iz (p_ch - p[ir,iz])/Rp0
    else # Bulk
        nϕ = ϕ[ir, iz+1]
        sϕ = ϕ[ir, iz-1]
        pp = p[ir, iz]
        pn = p[ir, iz+1]
        ps = p[ir, iz-1]
        dpz = (pn - ps) * 0.5dz1
        if sϕ <= 0 && nϕ <= 0 # North and south ghost cell: weird kink
            θz1 = ϕp/(ϕp - sϕ)
            θz2 = ϕp/(ϕp - nϕ)
            dpz = 0.5*dz1*(psub_l - pp)*(1/θz2 - 1/θz1)
        elseif sϕ <= 0 # South ghost cell
            θz = ϕp /(ϕp - sϕ)
            if θz >  θ_thresh
                dpz = (-psub_l/(1+θz)/θz + pp*(1-θz)/θz + pn*θz/(θz+1))*dz1 # Quadratic
                # dpz = (pp-psub_l)/θz*dz1 # Linear
            else
                dpz = (pn - psub_l)/(θz+1)*dz1 # Linear, further out
            end
        elseif nϕ <= 0 # North ghost cell
            θz = ϕp /(ϕp - nϕ)
            if θz >  θ_thresh
                dpz = ( psub_l/(θz+1)/θz - pp*(1-θz)/θz - ps*θz/(θz+1)) * dz1 # Quadratic extrapolation
                # dpz = (psub_l - pp)/θz*dz1 # Linear extrapolation
            else
                dpz = (psub_l - ps )/(θz+1)*dz1 #Linear, further out
            end
        else # No ghost cells
            dpz = (pn - ps) * 0.5*dz1 # Centered difference
        end

    end
    return dpr, dpz
end


"""
    compute_frontvel_mass(ϕ, T, p, dom::Domain, params; debug=false)

Generate an empty velocity field and compute velocity on `Γ⁺` (i.e. cells on Γ with ϕ>0). 
"""
function compute_frontvel_mass(u, T, p, dom::Domain, params; debug=false)

    @unpack k, ΔH, ρf, ϵ = params
    ϕ = ϕ_T_from_u(u, dom)[1]

    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
    b = eval_b(T, p, params)
     
    vf = zeros(eltype(ϕ), dom.nr, dom.nz, 2)
    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)


    for c in Γ⁺
        ir, iz = Tuple(c)
        dpr, dpz = compute_pderiv(u, T, p, ir, iz, dom, params)
        # if dpr > 0
        #     neighbors = [CI(i,j) for i in -1:1, j in 0]
        #     pnb = p[[c].+neighbors]
        #     bnb = b[[c].+neighbors]
        #     Tnb = T[[c].+neighbors]
        #     @warn "positive dpdr" c dpr pnb bnb Tnb
        # end

        # Boundary cases: use internal derivative
        if ir == dom.nr # Right boundary
            dϕdr = dϕdr_w[c]
        elseif ir == 1 # Left boundary
            dϕdr = dϕdr_e[c]
        else # Not at boundary: use according to pressure gradient
            dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        end
        if iz == dom.nz # Top boundary
            dϕdz = dϕdz_s[c]
        elseif iz == 1 # Bottom boundary
            dϕdz = dϕdz_n[c]
        else
            dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        end
        
        # # With extrapolation, no need for special boundary treatment
        # dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        # dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        # # Use dϕdr towards the interface, not according to pressure gradient
        # dϕdr = ((dϕdr_w[c] + dϕdr_e[c]) > 0 ? dϕdr_w[c] : dϕdr_e[c])
        # dϕdz = ((dϕdz_s[c] + dϕdz_n[c]) > 0 ? dϕdz_s[c] : dϕdz_n[c])
        
        # Normal is out of the ice
        # md = -b∇p⋅∇ϕ , >0 for sublimation occurring
        # v = md/ρ * -∇ϕ
        md_l = -b[c] * (dpr*dϕdr + dpz * dϕdz)
        vtot = md_l / ρf / ϵ 
        vf[c,1] = -vtot * dϕdr
        vf[c,2] = -vtot * dϕdz
        # if sign(vf[c,1]*dϕdr) != sign(vf[c,2]*dϕdz) && (dϕdr != 0 && dϕdz != 0)
        #     @warn "Velocity is messed up" vtot vf[c,1] dϕdr vf[c,2] dϕdz
        # end

        # if iz ∈ [2, 3]
        #     @info iz vf[c,2]
        # end

    end
    
    return vf, dϕdx_all
end

"""
    compute_frontvel_heat(u, T, p, dom::Domain, params; debug=false)

Generate an empty velocity field and compute velocity on `Γ⁺` (i.e. cells on Γ with ϕ>0). 
"""
function compute_frontvel_heat(u, T, dom::Domain, params; debug=false)

    @unpack k, ΔH, ρf, ϵ = params
    ϕ = ϕ_T_from_u(u, dom)[1]

    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
     
    vf = zeros(eltype(ϕ), dom.nr, dom.nz, 2)
    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    Qice = compute_Qice_nodry(u, T, dom, params)
    Q_ice_per_surf = Qice / compute_icesurf(ϕ, dom)

    for c in Γ⁺
        ir, iz = Tuple(c)
        dTr, dTz = compute_Tderiv(u, T, ir, iz, dom, params)

        # Boundary cases: use internal derivative
        if ir == dom.nr # Right boundary
            dϕdr = dϕdr_w[c]
        elseif ir == 1 # Left boundary
            dϕdr = dϕdr_e[c]
        else
            dϕdr = (dTr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        end
        if iz == dom.nz # Top boundary
            dϕdz = dϕdz_s[c]
        elseif iz == 1 # Bottom boundary
            dϕdz = dϕdz_n[c]
        else
            dϕdz = (dTz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        end

        # Normal is out of the ice
        # md =  , >0 for sublimation occurring
        # v = md/ρ * -∇ϕ
        q_grad = k* (dTr*dϕdr + dTz * dϕdz) 
        md_l = (q_grad + Q_ice_per_surf)/ΔH
        vtot = md_l / ρf / ϵ 
        vf[c,1] = -vtot * dϕdr
        vf[c,2] = -vtot * dϕdz

    end
    
    return vf, dϕdx_all
end

"""
    compute_frontvel_fixedspeed(v0, ϕ, dom::Domain)

Compute speed `v0` times normal vector on Γ⁺, positive side of interface.
"""
function compute_frontvel_fixedspeed(v0, ϕ, dom::Domain)
    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
     
    vf = zeros(eltype(ϕ), dom.nr, dom.nz, 2)
    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    for c in Γ⁺

        ir, iz = Tuple(c)
        # Boundary cases: use internal derivative
        if ir == dom.nr # Right boundary
            dϕdr = dϕdr_w[c]
        elseif ir == 1 # Left boundary
            dϕdr = dϕdr_e[c]
        else
            dϕdr = (dϕdr_w > 0 ? dϕdr_w[c] : dϕdr_e[c])
        end
        if iz == dom.nz # Top boundary
            dϕdz = dϕdz_s[c]
        elseif iz == 1 # Bottom boundary
            dϕdz = dϕdz_n[c]
        else
            dϕdz = (dϕdz_s > 0 ? dϕdz_s[c] : dϕdz_n[c])
        end

        # With extrapolation in WENO, no need for special boundary treatment
        # dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        # dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        
        vf[c,1] = -v0 * dϕdr
        vf[c,2] = -v0 * dϕdz

    end
    
    return vf, dϕdx_all
end
