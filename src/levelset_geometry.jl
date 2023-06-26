
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

function compute_icegl_area_weights(ϕ, dom)
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
    return zweights .* dom.dz
end

function compute_icesh_area_weights(ϕ, dom)
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
end


"""
    function compute_iceht_bottopcont(ϕ, dom)

Compute the height of ice at each point radially, as well as identify whether at system boundary or not.

Returns (`heights`, `bottom_contact`, `top_contact`), where 
- `heights` is a vector of floats
- `bottom_contact`, `top_contact` are vectors of bools
"""
function compute_iceht_bottopcont(ϕ, dom)
    bottom_contact = (ϕ[:,begin] .<= 0)
    top_contact = (ϕ[:,end] .<= 0)

    heights = fill(0.0, dom.nr)
    ice_contig = fill(true, dom.nr)
    for ir in axes(ϕ, 1)
        # Begin with 1 grid space per frozen cell
        heights[ir] += dom.dz * sum(ϕ[ir,:] .<= 0) 
        interfaces = 0
        for iz in 1:dom.nz-1
            ϕd, ϕu = ϕ[ir, iz:iz+1]
            if ϕd <= 0 && ϕu <= 0
                heights[ir] += dom.dz
            elseif ϕd <= 0
                θz = -(ϕd/(ϕd - ϕu))
                heights[ir] += dom.dz*θz
                interfaces += 1
            elseif ϕu <= 0
                θz = -(ϕu/(ϕu - ϕd))
                heights[ir] += dom.dz*θz
                interfaces += 1
            # else # Nothing needed to do in this case
            end
            if interfaces > 2 || (interfaces == 2 && ϕ[ir,begin] > 0)
                ice_contig[ir] = false
                @warn "Along a column, ice is not contiguous--more than two interfaces.
                    This case is not properly treated."
            end
        end
    end
    return heights, bottom_contact, top_contact
end

function compute_discrete_delta(i::Int, j::Int, ϕ, dom::Domain; ϵ=1e-10)
    rp_room = checkbounds(Bool, ϕ, i+1, j)
    rm_room = checkbounds(Bool, ϕ, i-1, j)
    zp_room = checkbounds(Bool, ϕ, i, j+1)
    zm_room = checkbounds(Bool, ϕ, i, j-1)
    if !rp_room
        D0r = (ϕ[i,j]-ϕ[i-1,j])*dom.dr1
    elseif !rm_room
        D0r = (ϕ[i+1,j]-ϕ[i,j])*dom.dr1
    else
        D0r = (ϕ[i+1,j]-ϕ[i-1,j])*0.5dom.dr1
    end
    if !zp_room
        D0z = (ϕ[i,j]-ϕ[i,j-1])*dom.dr1
    elseif !zm_room
        D0z = (ϕ[i,j+1]-ϕ[i,j])*dom.dr1
    else
        D0z = (ϕ[i,j+1]-ϕ[i,j-1])*0.5dom.dr1
    end

    norm∇ϕ = sqrt(D0r^2 + D0z^2 + ϵ)
    if rp_room && ϕ[i,j]*ϕ[i+1,j] <= 0
        D⁺r = (ϕ[i+1,j]-ϕ[i,j])*dom.dr1
        δr⁺ = abs(ϕ[i+1,j]*D0r)/abs(D⁺r)/norm∇ϕ*dom.dr1*dom.dz1
    else
        δr⁺ = 0
    end
    if rm_room && ϕ[i,j]*ϕ[i-1,j] < 0
        D⁻r = (ϕ[i,j]-ϕ[i-1,j])*dom.dr1
        δr⁻ = abs(ϕ[i-1,j]*D0r)/abs(D⁻r)/norm∇ϕ*dom.dr1*dom.dz1
    else
        δr⁻ = 0
    end
    if zp_room && ϕ[i,j+1]*ϕ[i,j] <= 0
        D⁺z = (ϕ[i,j+1]-ϕ[i,j])*dom.dr1
        δz⁺ = abs(ϕ[i,j+1]*D0z)/abs(D⁺z)/norm∇ϕ*dom.dr1*dom.dz1
    else
        δz⁺ = 0
    end
    if zm_room && ϕ[i,j]*ϕ[i,j-1] < 0
        D⁻z = (ϕ[i,j]-ϕ[i,j-1])*dom.dr1
        δz⁻ = abs(ϕ[i,j-1]*D0z)/abs(D⁻z)/norm∇ϕ*dom.dr1*dom.dz1
    else
        δz⁻ = 0
    end

    δ = δr⁺ + δr⁻ + δz⁺ + δz⁻
end

function compute_icesurf_δ(ϕ, dom)
    δ = compute_discrete_delta.(1:dom.nr, permutedims(1:dom.nz), [ϕ], [dom])
    SA = 2π*sum(δ .* dom.rgrid)*dom.dr*dom.dz
end