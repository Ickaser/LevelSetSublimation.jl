export compute_Qice, compute_icesurf
export compute_frontvel_withT, plot_frontvel

"""
    function compute_Qice(ϕ, dom::Domain, params)

Compute the total heat input into frozen domain from vial boundaries
"""
function compute_Qice(ϕ, dom::Domain, params)
    @unpack Q_sh, Q_ic, Q_gl = params
    # Heat flux from shelf, at bottom of vial
    botsurf = compute_icesh_area(ϕ, dom)
    Qbot = botsurf * Q_sh
    # Heat flux from glass, at outer radius
    outsurf = compute_icegl_area(ϕ, dom)
    Qout = outsurf * Q_gl
    # Volumetric heat throughout ice
    icevol = compute_icevol(ϕ, dom)
    Qvol = icevol * Q_ic
    return Qbot + Qout + Qvol
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
            θ = ϕpos / (ϕpos - ϕneg)
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
            dr = rmax - rmin # Absoluate value
            h = abs(zs[i+1] - zs[i])
            vol = 2π*h*(dr^2/3 + dr*rmin + rmax^2)

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
                    slope = (r2-r1)/(z2-z1)
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
            r2 = maximum(rs[i:i+1])
            r1 = minimum(rs[i:i+1])
            h = abs(zs[i+1] - zs[i])
            if h == 0 # Flat disk
                A = π * (r2^2 -r1^2)
            elseif r2 == 0 # Cylinder of radius 0
                A = 0
            elseif r1 == 0 # Cone outside
                l = hypot(r2, h)
                A = π*r2*l
            elseif r1 == r2
                A = π*r2*h
            else # Cone minus the top chunk
                h2 = h
                h1 = h2/r2/(1/r1 + 1/r2)
                h = h1 + h2
                l1 = hypot(h1, r1)
                l = hypot(h+h1, r2)
                A = π * (r2*l - r1*l1)
            end
            totsurf += A
        end
        rmids = [rs[i+1] + rs[i] for i in 1:length(rs)-1] .* 0.5
        dzs = abs.(zs[2:end] .- zs[1:end-1])
        # For surfaces which are exactly flat, replace dz with dr to get rdr
        # dzs[dzs .== 0] .= dr
        surf = sum(@. 2π*dzs*rmids) # Take a sum of conical surfaces
        totsurf += surf
    end
    return totsurf
end

"""

Quadratic ghost cell extrapolation (into frozen domain), second order finite differences, for T.
For ϕ derivatives, simple second order finite differences (one-sided at boundaries).
"""
function compute_heatflux(T, ϕ, ir::Int, iz::Int, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, nr, nz = dom
    @unpack k, Q_sh, Q_gl, Tf = params
    pT = T[ir, iz]
    pϕ = ϕ[ir, iz]
    
    if pϕ > 2dr || pϕ > 2dz || pϕ < -2dr || pϕ < -2dz
        @debug "Computing heat flux for cell which may not be at front." ir iz pϕ
    end

    # Enforce BCs explicitly for boundary cells
    if ir == 1 
        # dϕr = (ϕ[ir+1, iz] - pϕ) * dr1 # 1st order
        dϕr = (-0.5ϕ[ir+2, iz] +2ϕ[ir+1, iz] - 1.5pϕ) * dr1 # 2nd order 
        # dTr = (T[ir+1, iz] - pT) * dr1
        # dϕr = min(0, (ϕ[ir+1,iz] - pϕ)) * dr1 # Clamp to 0
        dTr = 0
    elseif ir == nr
        # dϕr = (pϕ - ϕ[ir-1, iz]) * dr1 # 1st order
        dϕr = (1.5pϕ - 2ϕ[ir-1, iz] + 0.5ϕ[ir-2,iz]) * dr1 # 2nd order
        # dTr = (pT - T[ir-1, iz]) * dr1
        # dϕr = min(0, (pϕ - ϕ[ir-1,iz])) * dr1 # Clamp to 0
        dTr = Q_gl / k
    else 
        # Bulk
        eϕ = ϕ[ir+1, iz]
        wϕ = ϕ[ir-1, iz]
        dϕr = (eϕ - wϕ) * 0.5dr1

        eT = T[ir+1, iz]
        wT = T[ir-1, iz]
        # West and east ghost cell: weird kink? Set to 0 and procrastinate
        if wϕ <= 0 && eϕ <= 0
            dϕr = (eϕ - wϕ) * 0.5*dr1 # Centered difference
            θr1 = pϕ/(pϕ - eϕ)
            θr2 = pϕ/(pϕ - wϕ)
            dTr = (Tf - pT)*(θr1 - θr2)/(θr2 * (2θr2-θr1))
        elseif wϕ <= 0 # West ghost cell
            θr = pϕ /(pϕ - wϕ)
            if θr > dr
                dTr = (-Tf/(1+θr)/θr + pT*(1-θr)/θr + eT*(θr)/(θr+1)) * dr1 # Quadratic extrapolation
            else 
                dTr = (eT - Tf)/(θr+1)*dr1 # Linear extrapolation from east
                # dTr = (eT - pT)       *dr1 # Linear extrapolation from east
                # @show pT-Tf θr+1
            end
        elseif eϕ <= 0 # East ghost cell
            θr = pϕ /(pϕ - eϕ)
            if θr > dr
                dTr = ( Tf/(θr+1)/θr - pT*(1-θr)/θr - wT*(3θr+1)/(θr+1)*0.25) * dr1 # Quadratic extrapolation
            else
                # eTg = (2Tf + (th-1)wT )/(th+1)
                dTr = (Tf - wT)/(θr+1)*dr1 # Linear extrapolation from west
            end
        else # No ghost cells
            dTr = (eT - wT) * 0.5*dr1 # Centered difference
        end
    end
            
    # Enforce BCs explicitly for boundary cells
    if iz == 1 
        # dϕz = (ϕ[ir, iz+1] - pϕ) * dr1 # 1st order
        dϕz = (-0.5ϕ[ir, iz+2] +2ϕ[ir, iz+1] - 1.5pϕ) * dz1 # 2nd order 
        dTz = Q_sh / k
    elseif iz == nz
        # dϕz = (pϕ - ϕ[ir, iz-1]) * dz1 # 1st order
        dϕz = (1.5pϕ - 2ϕ[ir, iz-1] + 0.5ϕ[ir,iz-2]) * dz1 # 2nd order
        # dϕz = min(0, (pϕ - ϕ[ir,iz-1])) * dz1 # Clamp to 0
        dTz = 0
    else 
        # Bulk
        nϕ = ϕ[ir, iz+1]
        sϕ = ϕ[ir, iz-1]
        dϕz = (nϕ - sϕ) * 0.5dz1

        nT = T[ir, iz+1]
        sT = T[ir, iz-1]
        # North and south ghost cell: weird kink
        if sϕ <= 0 && nϕ <= 0
            dϕz = (nϕ - sϕ) * 0.5*dz1 # Centered difference
            # dTz = 0
            θz1 = pϕ/(pϕ - nϕ)
            θz2 = pϕ/(pϕ - sϕ)
            dTz = (Tf - pT)*(θz1 - θz2)/(θz2 * (2θz2-θz1))
        elseif sϕ <= 0 # South ghost cell
            θz = pϕ /(pϕ - sϕ)
            if θz > dz
                dTz = (-Tf/(θz+1)/θz + pT*(1-θz)/θz + nT*θz/(θz+1)) * dz1 # Quadratic extrapolation
            else
                dTz = (nT - Tf)/(θz+1)*dz1
            end
        elseif nϕ <= 0 # North ghost cell
            θz = pϕ /(pϕ - nϕ)
            if θz > dz
                dTz = ( Tf/(θz+1)/θz - pT*(1-θz)/θz - sT*θz/(θz+1)) * dz1 # Quadratic extrapolation
            else
                dTz = (Tf - sT )/(θz+1)*dz1
            end
        else # No ghost cells
            dTz = (nT - sT) * 0.5*dz1 # Centered difference
        end
    end

    ngradϕ = hypot(dϕr, dϕz)
    dϕr /= ngradϕ
    dϕz /= ngradϕ
    return dϕr, dTr, dϕz, dTz
end

""" 
    compute_frontvel_withT(T, ϕ, ir::Int, iz::Int, dom::Domain, params, Qice_per_surf=nothing; debug=false)

Return `(vr, vz)` for the corresponding `(ir, iz)` location, based on temperature profile and z.

"""
function compute_frontvel_withT(T, ϕ, ir::Int, iz::Int, dom::Domain, params, Qice_per_surf=nothing; debug=false)
    # If Qice_surf (heat to ice divided by surface area) not supplied, compute it here
    if Qice_per_surf === nothing
        Qice = compute_Qice(ϕ, dom, params)
        icesurf = compute_icesurf(ϕ, dom)
        Qice_per_surf = Qice / icesurf
    end
    @unpack k, ΔH, ρf = params
    
    dϕr, dTr, dϕz, dTz = compute_heatflux(T, ϕ, ir, iz, dom, params)
    
    qtot = k*dTr*dϕr + k*dTz*dϕz + Qice_per_surf
    md = qtot / ΔH
    vtot = md / ρf
    
    return -vtot * dϕr, -vtot * dϕz

end

"""
    function plot_frontvel(ϕ, T, dom::Domain)

Calculate, then plot the front velocity given `ϕ` and `T`.

Meant for debugging, mostly. Scales all velocity arrows to have length 0.5.
Generates a freshplot().
"""
function plot_frontvel(ϕ, T, dom::Domain)
    front_cells = findall(identify_Γ(ϕ, dom) .& (ϕ .> 0))
    xs = []
    ys = []
    vrs = []
    vzs = []
    for cell in front_cells
        push!(xs, dom.rgrid[Tuple(cell)[1]])
        push!(ys, dom.zgrid[Tuple(cell)[2]])
        vr, vz = compute_frontvel_withT(T, ϕ, Tuple(cell)..., dom, T_params)
        push!(vrs, vr)
        push!(vzs, vz)
        # push!(vrs, get_front_vr(T, ϕ, Tuple(cell)..., T_params) )
        # push!(vzs, get_front_vz(T, ϕ, Tuple(cell)..., T_params) )
    end
    maxv = max(maximum(abs.(vrs)), maximum(abs.(vzs)))
    println("Maximum front velocity: $maxv")
    vrs ./= maxv * 2
    vzs ./= maxv * 2

    freshplot(dom)
    quiver!(xs, ys, quiver=(vrs, vzs))
end
# frontpos = findall(identify_Γ(fixed_ϕ) .& (fixed_ϕ .> 0))

