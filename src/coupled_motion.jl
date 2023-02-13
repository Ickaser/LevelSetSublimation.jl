export compute_Qice, compute_icesurf
export compute_frontvel_withT, plot_frontvel

"""
    function compute_Qice(ϕ, dom::Domain, params)

Compute the total heat input into frozen domain from vial boundaries

This function is currently not treating the sharp interface carefully.
TODO
"""
function compute_Qice(ϕ, dom::Domain, params)
    @unpack Q_sh, Q_ic, Q_gl = params

    # Heat flux from shelf, at bottom of vial
    botϕ = ϕ[:,1]
    botr = dom.rgrid[botϕ .<=0]
    botsurf = 2π * sum(botr) * dom.dr # Should interpolate at edge.
    Qbot = botsurf * Q_sh

    # Heat flux from glass, at outer radius
    outϕ = ϕ[end,:]
    outsurf = 2π * dom.rmax * dom.dz * sum(outϕ .<=0)
    Qout = outsurf * Q_gl

    # Volumetric heat throughout ice: compute like bottom surface above
    icevol = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π
    Qvol = icevol * Q_ic
    
    # println("Qbot = $Qbot, Qout = $Qout")
    return Qbot + Qout + Qvol
end


"Geometric: compute area of cone-shaped sections of interface."
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
    compute_frontvel_withT(T, ϕ, ir, iz, dom::Domain, params, Qice_per_surf=nothing; debug=false)

Return (vr, vz) for the corresponding (ir, iz) location

We must take temperature derivatives into positive ϕ, since in negative ϕ (frozen) we have constant T.
This means that temperature is computed "downwind" in the level-set-reinitialization sense.
To get the normal vector, take "upwind" differences (but assume positive ϕ).
(Since this isn't the RHS of a hyperbolic PDE, it's less important to avoid downwind?)
First order (downwind/upwind) differences, currently, to simplify reasoning.
"""
function compute_frontvel_withT(T, ϕ, ir, iz, dom::Domain, params, Qice_per_surf=nothing; debug=false)
    # dr = dom.dr
    # dz = dom.dz
    # nr = dom.nr
    # nz = dom.nz
    @unpack dr, dz, dr1, dz1, nr, nz = dom
        
    # If Qice_surf (heat to ice divided by surface area) not supplied, compute it from shelf
    # This should be outside this function. TODO
    if Qice_per_surf === nothing
        Qice = compute_Qice(ϕ, dom, params)
        icesurf = compute_icesurf(ϕ, dom)
        Qice_per_surf = Qice / icesurf
    end
    
    @unpack k, ΔH, ρf, Q_sh, Q_gl = params
    pT = T[ir, iz]
    pϕ = ϕ[ir, iz]
    
    if pϕ > 2dr || pϕ > 2dz || pϕ < -2dr || pϕ < -2dz
        @warn "Computing front velocity for cell which may not be at front."
    end

    # Enforce BCs explicitly for boundary cells
    if ir == 1 
        dϕr = (ϕ[ir+1, iz] - pϕ) * dr1 # Possibly downwind
        # dTr = (T[ir+1, iz] - pT) * dr1
        # dϕr = min(0, (ϕ[ir+1,iz] - pϕ)) * dr1 # Clamp to 0
        dTr = 0
    elseif ir == nr
        dϕr = (pϕ - ϕ[ir-1, iz]) * dr1 # Possibly downwind
        # dTr = (pT - T[ir-1, iz]) * dr1
        # dϕr = min(0, (pϕ - ϕ[ir-1,iz])) * dr1 # Clamp to 0
        dTr = Q_gl / k
    else 
        # In the bulk, take difference towards positive ϕ
        eϕ = ϕ[ir+1, iz]
        wϕ = ϕ[ir-1, iz]
        if eϕ > pϕ && wϕ > pϕ
            dϕr = (eϕ - wϕ) * 0.5*dr1 # Centered difference
            dTr = (T[ir+1,iz] - T[ir-1,iz]) * 0.5*dr1 # Centered difference
        elseif eϕ > pϕ
            dϕr = (pϕ         - wϕ)*dr1
            dTr = (T[ir+1,iz] - pT)*dr1
        elseif wϕ > pϕ
            dϕr = (eϕ -         pϕ)*dr1
            dTr = (pT - T[ir-1,iz])*dr1
        else
            dϕr = (eϕ - wϕ) * 0.5*dr1 # Centered difference
            dTr = (T[ir+1,iz] - T[ir-1,iz]) * 0.5*dr1 # Centered difference
            # dϕr = 0
            # dTr = 0 # Technically would take a central difference, but is irrelevant
        end
    end
            
    # Enforce BCs explicitly for boundary cells
    if iz == 1 
        # dϕz = (ϕ[ir, iz+1] - pϕ) * dz1
        # dTz = (T[ir, iz+1] - pT) * dz1
        dϕz = min(0, (ϕ[ir, iz+1] - pϕ)) * dz1 # Clamp to 0 if downwind
        dTz = - Q_sh / k
    elseif iz == nz
        # dϕz = (pϕ - ϕ[ir, iz-1]) * dz1
        # dTz = (pT - T[ir, iz-1]) * dz1
        dϕz = max(0, (pϕ - ϕ[ir, iz-1])) * dz1 # Clamp to 0 if downwind
        dTz = 0
    else 
        # In the bulk, take difference towards positive ϕ
        nϕ = ϕ[ir, iz+1]
        sϕ = ϕ[ir, iz-1]
        if nϕ > pϕ && sϕ > pϕ
            dϕz = (nϕ         - sϕ)*dz1 # Centered
            dTz = (T[ir,iz+1] - T[ir,iz-1])*0.5*dz1 # Centered
        elseif nϕ > pϕ
            dϕz = (pϕ -         sϕ)*dz1 # Downward
            dTz = (T[ir,iz+1] - pT)*dz1 # Upward
        elseif sϕ > pϕ
            dϕz = (nϕ -         pϕ)*dz1 # Upward
            dTz = (pT - T[ir,iz-1])*dz1 # Downward
        else
            # dϕz = 0
            # dTz = 0 # Technically would take a central difference, but is irrelevant
            dϕz = (nϕ - sϕ)*0.5*dz1 # Centered
            dTz = (T[ir,iz+1] - T[ir,iz-1])*0.5*dz1 # Centered
        end
    end
    
    # q = k*(abs(dTr*dϕr) + abs(dTz*dϕz)) # Assumed: heat is going into frozen domain.
    # q += Qice_per_surf
    # md = q/ΔH
    # vtot = md / ρf

    # qr = k*dTr + Qice_per_surf * dϕr
    # qz = k*dTz + Qice_per_surf * dϕz

    if debug
        println("dTr = $dTr, dTz = $dTz, Qice = $Qice_per_surf, c=$((ir,iz))")
    end

    # qr = -k* dTr , qtot = -q⋅n = -qr*dϕr - qz*dϕz
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

