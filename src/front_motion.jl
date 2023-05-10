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

    # @info "Q" Q_glshvol Qgl Qsub

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
    @unpack Kv, Kgl, Q_ic, Q_ck, Tsh= params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)

    # Heat flux from shelf, at bottom of vial

    # Approximately compute area
    # rweights = zeros(Float64, dom.nr) 
    # for ir in 1:dom.nr-1
    #     ϕl, ϕr = ϕ[ir:ir+1, 1]
    #     if ϕl <= 0 && ϕr <= 0
    #         rweights[ir:ir+1] .+= 0.5
    #     elseif ϕl <= 0
    #         rweights[ir] += -(ϕl/(ϕl - ϕr))*0.5
    #     elseif ϕr <= 0
    #         rweights[ir+1] += -(ϕr/(ϕr - ϕl))*0.5
    #     # else # Nothing needed to do in this case
    #     end
    # end
    # Qsh = 2π* Kv * dom.dr* sum(rweights .* dom.rgrid .* (Tsh .- T[:,1] ))

    # Explicitly compute area
    rweights = zeros(Float64, dom.nr) 
    for ir in 1:dom.nr-1
        ϕl, ϕr = ϕ[ir:ir+1, 1]
        if ϕl <= 0 && ϕr <= 0
            rmid = (dom.rgrid[ir] + dom.rgrid[ir+1])/2
            rweights[ir] += rmid^2/2
            rweights[ir+1] -= rmid^2/2
        elseif ϕl <= 0
            θr = -(ϕl/(ϕl - ϕr))
            rmid =dom.rgrid[ir] + dom.dr*θr
            rweights[ir] += rmid^2/2
        elseif ϕr <= 0
            θr = -(ϕr/(ϕr - ϕl))
            rmid = dom.rgrid[ir+1] - dom.dr*θr
            rweights[ir+1] -= rmid^2/2
        # else # Nothing needed to do in this case
        end
    end
    if ϕ[dom.nr] <= 0
        rweights[dom.nr] += dom.rmax^2/2
    end

    Qsh = 2π* Kv * sum(rweights .* (Tsh .- T[:,1] ))
    

    # Heat flux from glass, at outer radius
    zweights = zeros(Float64, dom.nz) 
    for iz in 1:dom.nz-1
        ϕl, ϕr = ϕ[end, iz:iz+1]
        if ϕl <= 0 && ϕr <= 0
            zweights[iz:iz+1] .+= 0.5
        elseif ϕl <= 0
            zweights[iz] += -(ϕl/(ϕl - ϕr))*0.5
        elseif ϕr <= 0
            zweights[iz+1] += -(ϕr/(ϕr - ϕl))*0.5
        # else # Nothing needed to do in this case
        end
    end
    Qgl = 2π*dom.rmax * Kgl * dom.dz * sum(zweights .* ( Tgl .- T[end,:]))

    # Volumetric heat in cake and ice
    icevol = compute_icevol(ϕ, dom)
    Qvol = icevol * Q_ic 

    # @info "geom" icevol get_subf_r(ϕ, dom)
    return Qsh + Qgl + Qvol
end

"""
    compute_topmassflux(ϕ, T, p, dom::Domain, params)

Compute total mass flow through top of the cake (that is, mass flux integrated across top surface).
"""
function compute_topmassflux(u, T, p, dom::Domain, params)
    dpdz = [compute_pderiv(u, T, p, ir, dom.nz, dom, params)[2] for ir in 1:dom.nr]
    b = eval_b(T, p, params)[:,end] # all r, top of z
    md = - 2π * sum(dpdz .*b .* dom.rgrid) * dom.dr
    # @info "massflux" dpdz b md 1/dom.dz
    return md
end

"""
    compute_icesh_area(ϕ, dom::Domain)

Compute area of where ice meets bottom surface of vial.
"""
function compute_icesh_area(ϕ, dom::Domain)
    botϕ = ϕ[:,1]
    # botr = dom.rgrid[botϕ .> 0]
    # botsurf = 2π * sum(botr) * dom.dr # Not interpolating at edge.
    botsurf::Float64 = 0.0
    for i in 1:length(botϕ)-1
        ϕl, ϕr = botϕ[i:i+1]
        if ϕl <= 0 && ϕr <= 0 # full segment in ice
            rmin = dom.rgrid[i]
            rmax = dom.rgrid[i+1]
        elseif ϕl > 0 && ϕr <= 0 # right of segment in ice
            θ = ϕr / (ϕr - ϕl)
            rmin = dom.rgrid[i] + θ*dom.dr
            rmax = dom.rgrid[i+1]
        elseif ϕl <= 0 && ϕr > 0 # left of segment in ice
            θ = ϕl / (ϕl - ϕr)
            rmin = dom.rgrid[i] 
            rmax = dom.rgrid[i+1] + θ*dom.dr
        else
            continue # Full segment is dried
        end
        botsurf += π*(rmax^2 - rmin^2)
    end
    return botsurf
end

"""
    compute_icegl_area(ϕ, dom::Domain)

Compute area of where ice meets radial outer surface of vial.
"""
function compute_icegl_area(ϕ, dom::Domain)
    outϕ = ϕ[end,:]
    # outsurf = dom.rmax*dom.dz*2π * sum(outϕ .> 0)
    outsurf::Float64 = 0.0
    outcoeff = dom.rmax*dom.dz*2π
    for seg in zip(outϕ[begin:end-1], outϕ[begin+1:end])
        if seg[1] <= 0 && seg[2] <= 0 # full segment in ice
            outsurf += outcoeff
        elseif (seg[1] <= 0) ⊻ (seg[2] <= 0)  # ⊻ xor: part of segment in ice
            ϕneg, ϕpos = minmax(seg...)
            θ = ϕneg / (ϕneg - ϕpos) # Inverted from our ghost cells- we want the negative part
            outsurf += θ*outcoeff
        end
    end
    return outsurf
end

"""
    compute_icevol(ϕ, dom::Domain)

Compute volume of ice (region where ϕ < 0) as sum of cylinders and cones.

Cylinders at outer edge (maximum radius), cones everywhere else, which are identified by interface.
"""
function compute_icevol(ϕ, dom::Domain)
    # One-line version, without any cones:
    # icevol = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π

    function loc_in_grid(r, z, dom)
        rloc = findlast(r .>= dom.rgrid)
        zloc = findlast(z .>= dom.zgrid)
        return CI(rloc, zloc)
    end

    totvol::Float64 = 0.0
    # Volume from outer cylinder edge: R/2 * outer surface area
    totvol += dom.rmax/2 * compute_icegl_area(ϕ, dom)

    cl = contour(dom.rgrid,dom.zgrid,ϕ, 0)
    for line in lines(cl) # iterate over disconnected segments of contour
        rs, zs = coordinates(line) # coordinates of this section of contour 
        for i in 1:length(rs)-1

            r1, r2 = rs[i:i+1]
            z1, z2 = zs[i:i+1]
            # Compute volume of cone slice: cone minus the top chunk
            rmin, rmax = minmax(r1, r2)
            dr = rmax - rmin # Absolute value
            if dr != 0 # Cone
                h = abs(z2 - z1)
                H = rmax*h/dr
                vol = π/3 * (H*rmax^2 - (H-h)*rmin^2)
            else # Cylinder
                vol = π*rmax^2*abs(z2-z1)
            end
            # If cone contains ice, add to volume.
            # If cone contains air, subtract, because surrounded by other ice.
            # Six cases of possible edge connections in a box to consider; worked out in my tablet notes
            ij1 = loc_in_grid(r1, z1, dom)
            ij2 = loc_in_grid(r2, z2, dom)
            mi = min(ij1, ij2) # Bottom left index
            m = max(ij1, ij2) - mi # Local max of indices: (0,0), (0,1), (1,0), or (1,1)
            if m == CI(1, 1) || m == CI(0, 0) # Case 1 or 6 from notes: top edge to right edge, left edge to bottom edge
                if ϕ[mi] <= 0
                    totvol += vol
                else
                    totvol -= vol
                end
            elseif m == CI(0, 1) # Case 2 and 3 from notes, handled together: top to bottom or top to left
                if ϕ[mi+m] < 0
                    totvol += vol
                else
                    totvol -= vol
                end
            elseif m == CI(1, 0) # Case 4 and 5 from notes, separate
                if ϕ[mi+CI(0,1)] * ϕ[mi] > 0 # Case 4: right to bottom
                    if ϕ[mi] < 0
                        totvol += vol
                    else
                        totvol -= vol
                    end
                else   # Case 5: left to right, slope matters for in or out
                    slope = (z2-z1)/(r2-r1)
                    totvol += vol * sign(ϕ[mi] * slope)
                end
            else
                @warn "Didn't get a case match for volume" ij1 ij2 vol
            end
        end
    end
    return totvol
end

"""
    compute_icesurf(ϕ, dom::Domain)

Compute the surface area of interface, treating as cone-shaped sections.
"""
function compute_icesurf(ϕ, dom::Domain)
    totsurf = 0.0
    # cl = levels(contours(dom.rgrid,dom.zgrid,ϕ, 0))[1]
    cl = contour(dom.rgrid,dom.zgrid,ϕ, 0)
    for line in lines(cl)
        rs, zs = coordinates(line) # coordinates of this line segment
        for i in 1:length(rs)-1
            # Compute area: cone minus the top chunk
            r1, r2 = extrema(rs[i:i+1])
            h = abs(zs[i+1] - zs[i])
            if h == 0 # Flat disk
                A = π * (r2^2 -r1^2)
            elseif r2 == 0 # Cylinder of radius 0
                A = 0
            elseif r1 == 0 # Cone outside
                l = hypot(r2, h)
                A = π*r2*l
            elseif (r2 - r1)/(h) < 1e-10 # Cylincrical, or at least nearly. Cone formula is singular here
                A = 2*π*r2*h
            else # Cone minus the top chunk
                # 1: top chunk
                # 2: bottom chunk (part of interest)
                l2 = hypot(r2-r1, h)
                l1 = l2/r2 / (1/r1 - 1/r2)
                l = l1 + l2
                A = π * (r2*(l1+l2) - r1*l1)
            end
            totsurf += A
            # @show A
        end
    end
    return totsurf
end

"""
    compute_Tderiv(u, T, ir::Int, iz::Int, dom::Domain, params)

Compute (`∂T/∂r`, `∂T/∂z`) at point `(ir, iz)`, given system state `u`, `T`.

Quadratic ghost cell extrapolation (into frozen domain), second order finite differences, for T.
"""
function compute_Tderiv(u, T, ir::Int, iz::Int, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, nr, nz = dom
    @unpack k, Kv, Kgl, Tsh = params

    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    pT = T[ir, iz]
    ϕp = ϕ[ir, iz]
    r = dom.rgrid[ir]
    
    if ϕp > 2dr || ϕp > 2dz || ϕp < -2dr || ϕp < -2dz
        @debug "Computing heat flux for cell which may not be at front." ir iz ϕp
    end

    θ_thresh = max(1/nr, 1/nz)

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
        # West and east ghost cell: weird kink? Set to 0 and procrastinate
        if wϕ <= 0 && eϕ <= 0
            θr1 = ϕp/(ϕp - wϕ)
            θr2 = ϕp/(ϕp - eϕ)
            dTr = 0.5*dr1*(Tf - pT)*(1/θr2 - 1/θr1)
        elseif wϕ <= 0 # West ghost cell
            θr = ϕp /(ϕp - wϕ)
            if θr > θ_thresh
                # dTr = (-Tf/(1+θr)/θr + pT*(1-θr)/θr + eT*(θr)/(θr+1)) * dr1 # Quadratic extrapolation
                dTr = (pT - Tf)/θr * dr1 # LInear extrapolation
            else 
                dTr = (eT - Tf)/(θr+1)*dr1 # Linear extrapolation from east
            end
        elseif eϕ <= 0 # East ghost cell
            θr = ϕp /(ϕp - eϕ)
            if θr > θ_thresh
                # dTr = ( Tf/(θr+1)/θr - pT*(1-θr)/θr - wT*θr/(θr+1)) * dr1 # Quadratic extrapolation
                dTr = (Tf - pT)/θr * dr1 # LInear extrapolation
            else
                dTr = (Tf - wT)/(θr+1)*dr1 # Linear extrapolation from west
            end
        else # No ghost cells
            dTr = (eT - wT) * 0.5*dr1 # Centered difference
        end
    end
            
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
            dTz = 0.5*dz1*(Tf - pT)*(1/θz2 - 1/θz1)
        elseif sϕ <= 0 # South ghost cell
            θz = ϕp /(ϕp - sϕ)
            if θz > θ_thresh
                # dTz = (-Tf/(θz+1)/θz + pT*(1-θz)/θz + nT*θz/(θz+1)) * dz1 # Quadratic extrapolation
                dTz = (pT - Tf)/θz * dz1 # Linear extrapolation
            else
                dTz = (nT - Tf)/(θz+1)*dz1
            end
        elseif nϕ <= 0 # North ghost cell
            θz = ϕp /(ϕp - nϕ)
            if θz > θ_thresh
                # dTz = ( Tf/(θz+1)/θz - pT*(1-θz)/θz - sT*θz/(θz+1)) * dz1 # Quadratic extrapolation
                dTz = (Tf - pT)/θz * dz1 # Linear extrapolation
            else
                dTz = (Tf - sT )/(θz+1)*dz1
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
function compute_pderiv(u, T, p, ir::Int, iz::Int, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, nr, nz = dom
    @unpack p_ch, Rp0 = params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    p_sub = calc_psub(Tf)
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
            θr1 = ϕp/(ϕp - wϕ)
            θr2 = ϕp/(ϕp - eϕ)
            dpr = 0.5*dr1*(p_sub - pp)*(1/θr2 - 1/θr1)
        elseif wϕ <= 0 # West ghost cell
            θr = ϕp /(ϕp - wϕ)
            if θr > θ_thresh
                dpr = (-p_sub/(1+θr)/θr + pp*(1-θr)/θr + pe*(θr)/(θr+1)) * dr1 # Quadratic extrapolation
                # dpr = (pp - p_sub)/θr*dr1 # Linear extrapolation
            else 
                dpr = (pe - p_sub)/(θr+1)*dr1 # Linear extrapolation, further out
            end
        elseif eϕ <= 0 # East ghost cell
            θr = ϕp /(ϕp - eϕ)
            if θr > θ_thresh
                dpr = ( p_sub/(θr+1)/θr - pp*(1-θr)/θr - pw*θr/(θr+1)) * dr1 # Quadratic extrapolation
                # dpr = (p_sub - pp)/θr*dr1 # Linear extrapolation
            else
                dpr = (p_sub - pw)/(θr+1)*dr1 # Linear extrapolation, further out
            end
        else # No ghost cells
            dpr = (pe - pw) * 0.5*dr1 # Centered difference
        end
        # if (ir, iz) == (43, 38)
        # @info "pderiv" ir iz wϕ ϕp eϕ  pw pp pe dpr
        # end

    end
            
    if iz == 1 
        # Enforce BCs explicitly for Neumann boundary cells
        dpz = 0
    elseif iz == nz 
        # Robin boundary condition: employ explicitly
        b = eval_b(T[ir,iz], p[ir,iz], params)
        if length(b) == 1
            bp = b
        else
            bp = b[ir,iz]
        end
        dpz = Rp0/bp*(p_ch - p[ir,iz]) 
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
            dpz = 0.5*dz1*(p_sub - pp)*(1/θz2 - 1/θz1)
        elseif sϕ <= 0 # South ghost cell
            θz = ϕp /(ϕp - sϕ)
            if θz >  θ_thresh
                dpz = ((-p_sub + pp*(1-θz))/θz + pn)*0.5*dz1 # Quadratic
                # dpz = (pp-p_sub)/θz*dz1 # Linear
            else
                dpz = (pn - p_sub)/(θz+1)*dz1 # Linear, further out
            end
        elseif nϕ <= 0 # North ghost cell
            θz = ϕp /(ϕp - nϕ)
            if θz >  θ_thresh
                dpz = ( p_sub/(θz+1)/θz - pp*(1-θz)/θz - ps*θz/(θz+1)) * dz1 # Quadratic extrapolation
                # dpz = (p_sub - pp)/θz*dz1 # Linear extrapolation
            else
                dpz = (p_sub - ps )/(θz+1)*dz1 #Linear, further out
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
    dϕdr_e, dϕdr_w, dϕdz_n, dϕdz_s = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    for c in Γ⁺
        ir, iz = Tuple(c)
        dpr, dpz = compute_pderiv(u, T, p, ir, iz, dom, params)

        # # Boundary cases: use internal derivative
        # if ir == dom.nr # Right boundary
        #     dϕdr = dϕdr_w[c]
        # elseif ir == 1 # Left boundary
        #     dϕdr = dϕdr_e[c]
        # else
        #     dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        # end
        # if iz == dom.nz # Top boundary
        #     dϕdz = dϕdz_s[c]
        # elseif iz == 1 # Bottom boundary
        #     dϕdz = dϕdz_n[c]
        # else
        #     dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        # end
        
        # With extrapolation, no need for special boundary treatment
        dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        
        # Normal is out of the ice
        # md = -b∇p⋅∇ϕ , >0 for sublimation occurring
        # v = md/ρ * -∇ϕ
        md_l = -b[c] * (dpr*dϕdr + dpz * dϕdz)
        vtot = md_l / ρf / ϵ 
        vf[c,1] = -vtot * dϕdr
        vf[c,2] = -vtot * dϕdz

    end
    
    return vf, dϕdx_all
end

"""
    compute_frontvel_heat(ϕ, T, p, dom::Domain, params; debug=false)

Generate an empty velocity field and compute velocity on `Γ⁺` (i.e. cells on Γ with ϕ>0). 
"""
function compute_frontvel_heat(u, T, dom::Domain, params; debug=false)

    @unpack k, ΔH, ρf, ϵ = params
    ϕ = ϕ_T_from_u(u, dom)[1]

    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
     
    vf = zeros(eltype(ϕ), dom.nr, dom.nz, 2)
    dϕdr_e, dϕdr_w, dϕdz_n, dϕdz_s = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    Qice = compute_Qice_nodry(u, T, dom, params)
    Q_ice_per_surf = Qice / compute_icesurf(ϕ, dom)

    for c in Γ⁺
        ir, iz = Tuple(c)
        dTr, dTz = compute_Tderiv(u, T, ir, iz, dom, params)

        # # Boundary cases: use internal derivative
        # if ir == dom.nr # Right boundary
        #     dϕdr = dϕdr_w[c]
        # elseif ir == 1 # Left boundary
        #     dϕdr = dϕdr_e[c]
        # else
        #     dϕdr = (dTr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        # end
        # if iz == dom.nz # Top boundary
        #     dϕdz = dϕdz_s[c]
        # elseif iz == 1 # Bottom boundary
        #     dϕdz = dϕdz_n[c]
        # else
        #     dϕdz = (dTz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        # end

        # With extrapolation, no need for special boundary treatment
        dϕdr = (dTr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        dϕdz = (dTz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        
        # Normal is out of the ice
        # md =  , >0 for sublimation occurring
        # v = md/ρ * -∇ϕ
        q_grad = k* (dTr*dϕdr + dTz * dϕdz) 
        # @info "aha"  q_grad Qice Q_ice_per_surf
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
    dϕdr_e, dϕdr_w, dϕdz_n, dϕdz_s = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    for c in Γ⁺

        # ir, iz = Tuple(c)
        # # Boundary cases: use internal derivative
        # if ir == dom.nr # Right boundary
        #     dϕdr = dϕdr_w[c]
        # elseif ir == 1 # Left boundary
        #     dϕdr = dϕdr_e[c]
        # else
        #     dϕdr = (dϕdr_w > 0 ? dϕdr_w[c] : dϕdr_e[c])
        # end
        # if iz == dom.nz # Top boundary
        #     dϕdz = dϕdz_s[c]
        # elseif iz == 1 # Bottom boundary
        #     dϕdz = dϕdz_n[c]
        # else
        #     dϕdz = (dϕdz_s > 0 ? dϕdz_s[c] : dϕdz_n[c])
        # end

        # With extrapolation in WENO, no need for special boundary treatment
        dϕdr = (dpr < 0 ? dϕdr_w[c] : dϕdr_e[c])
        dϕdz = (dpz < 0 ? dϕdz_s[c] : dϕdz_n[c])
        
        vf[c,1] = -v0 * dϕdr
        vf[c,2] = -v0 * dϕdz

    end
    
    return vf, dϕdx_all
end
