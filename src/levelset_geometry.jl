
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
    return rweights
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
        # # Begin with 1 grid space per frozen cell
        # heights[ir] += dom.dz * sum(ϕ[ir,:] .<= 0) 
        interfaces = 0
        for iz in 1:dom.nz-1
            ϕd, ϕu = ϕ[ir, iz:iz+1]
            if ϕd <= 0 && ϕu <= 0
                heights[ir] += dom.dz
            elseif ϕd <= 0
                θz = (ϕd/(ϕd - ϕu))
                heights[ir] += dom.dz*θz
                interfaces += 1
            elseif ϕu <= 0
                θz = (ϕu/(ϕu - ϕd))
                heights[ir] += dom.dz*θz
                interfaces += 1
            # else # Nothing needed to do in this case
            end
        end
        if interfaces > 2 || (interfaces == 2 && bottom_contact[ir])
            ice_contig[ir] = false
            @warn "Along a column, ice is not contiguous--more than two interfaces.
                This case is not properly treated." interfaces 
            println(ϕ[ir, :] .<= 0)
                
        end
    end
    if !any(bottom_contact)
        display(heat(ϕ, dom))
    end
    return heights, bottom_contact, top_contact
end

function compute_icesurf_δ(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    SA = 2π*sum(δ .* dom.rgrid)*dom.dr*dom.dz
end
function compute_icevol_H(ϕ, dom)
    H = compute_discrete_H(ϕ, dom)
    SA = 2π*sum(H .* dom.rgrid)*dom.dr*dom.dz
end

# --------------- Min & Gibou 2008 on surface integrals

Pij(Pi, Pj, ϕi, ϕj) = Pi*ϕj/(ϕj-ϕi) + Pj*ϕi/(ϕi-ϕj)
Pi(i, j, dom) = [dom.rgrid[i], dom.zgrid[j]]
Pi(ij, dom) = [dom.rgrid[Tuple(ij)[1]], dom.zgrid[Tuple(ij)[2]]]
distance(Pi, Pj) = norm(Pi.-Pj, 2)
area(Pi, Pj, Pk) = abs(Pi[1]*(Pj[2]-Pk[2]) + Pj[1]*(Pk[2]-Pi[2]) + Pk[1]*(Pi[2]-Pj[2]))/2

"""
    calc_δ0(ϕ0, ϕ1, ϕ2, P0, P1, P2)

Implementation of Table 1 from Min & Gibou, 2008.
"""
function calc_δ0(ϕ0, ϕ1, ϕ2, P0, P1, P2)
    s0, s1, s2 = sign.([ϕ0, ϕ1, ϕ2])
    if s0 == s1 && s1 == s2
        δ0 = 0
    elseif s0 == s1 && s0 != s2
        P02 = Pij(P0, P2, ϕ0, ϕ2)
        P12 = Pij(P1, P2, ϕ1, ϕ2)
        δ0 = distance(P02, P12)/2 * ϕ2/(ϕ2-ϕ0)
    elseif s0 != s1 && s0 == s2
        P01 = Pij(P0, P1, ϕ0, ϕ1)
        P21 = Pij(P2, P1, ϕ2, ϕ1)
        δ0 = distance(P01, P21)/2*ϕ1/(ϕ1-ϕ0)
    elseif s0 != s1 && s0 != s2
        P01 = Pij(P0, P1, ϕ0, ϕ1)
        P02 = Pij(P0, P2, ϕ0, ϕ2)
        δ0 = distance(P01, P02)/2 * (ϕ1/(ϕ1-ϕ0) + ϕ2/(ϕ2-ϕ0))
    else
        @warn "Uncaught case somehow"
    end
    return δ0
end

"""
    compute_local_δ(ij::CartesianIndex, ϕ, dom)

Return the discrete delta for level set `ϕ` at location `ij`.

Implementation of the final expression for a discrete delta in Min & Gibou, 2008.
"""
function compute_local_δ(ij::CartesianIndex, ϕ, dom)
    δij = 0
    simplex_vec = []
    if checkbounds(Bool, ϕ, ij + CI(1,1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(1,0), CI(1,1)])
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(0,1), CI(1,1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(-1, -1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(-1,0), CI(-1,-1)])
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(0,-1), CI(-1,-1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(1, -1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(1,0), CI(0,-1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(-1, 1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(-1,0), CI(0,1)])
    end
    for simplex in simplex_vec
        δij += calc_δ0(ϕ[simplex]..., Pi.(simplex, [dom])...)
    end
    return δij*dom.dr1*dom.dz1
end

"""
    compute_discrete_δ(ϕ, dom)

Compute the discrete Dirac δ across the domain, for use in surface integrals. 
"""
function compute_discrete_δ(ϕ, dom)
    compute_local_δ.(CartesianIndices(ϕ), [ϕ], [dom])
end

"""
    calc_H0(ϕ0, ϕ1, ϕ2, P0, P1, P2)

Implementation of Table 1 from Min & Gibou, 2008.
"""
function calc_H0(ϕ0, ϕ1, ϕ2, P0, P1, P2)
    s0, s1, s2 = sign.([ϕ0, ϕ1, ϕ2])
    if s0 == 1
        return area(P0, P1, P2)/3 - calc_H0(-ϕ0, -ϕ1, -ϕ2, P0, P1, P2)
    end

    if s0 == s1 && s1 == s2
        H0 = area(P0, P1, P2)/3
    elseif s0 == s1 && s0 != s2
        P02 = Pij(P0, P2, ϕ0, ϕ2)
        P12 = Pij(P1, P2, ϕ1, ϕ2)
        H0 = area(P0, P1, P2)/3 - area(P02, P12, P2)/3 * ϕ2/(ϕ2-ϕ0)
    elseif s0 != s1 && s0 == s2
        P01 = Pij(P0, P1, ϕ0, ϕ1)
        P21 = Pij(P2, P1, ϕ2, ϕ1)
        H0 = area(P0, P1, P2)/3 - area(P01, P21, P1)/3 * ϕ1/(ϕ1-ϕ0)
    elseif s0 != s1 && s0 != s2
        P01 = Pij(P0, P1, ϕ0, ϕ1)
        P02 = Pij(P0, P2, ϕ0, ϕ2)
        H0 = area(P01, P02, P0)/3 * (1 + ϕ1/(ϕ1-ϕ0) + ϕ2/(ϕ2-ϕ0))
    else
        @warn "Uncaught case somehow"
    end
    return H0
end
"""
    compute_local_H(ij::CartesianIndex, ϕ, dom)

Return the discrete Heaviside for level set `ϕ` at location `ij`.

Implementation of the final expression for a discrete delta in Min & Gibou, 2008.
"""
function compute_local_H(ij::CartesianIndex, ϕ, dom)
    Hij = 0
    simplex_vec = []
    if checkbounds(Bool, ϕ, ij + CI(1,1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(1,0), CI(1,1)])
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(0,1), CI(1,1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(-1, -1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(-1,0), CI(-1,-1)])
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(0,-1), CI(-1,-1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(1, -1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(1,0), CI(0,-1)])
    end
    if checkbounds(Bool, ϕ, ij + CI(-1, 1))
        push!(simplex_vec, [ij] .+ [CI(0,0), CI(-1,0), CI(0,1)])
    end
    for simplex in simplex_vec
        # @info "simplex" simplex calc_H0(ϕ[simplex]..., Pi.(simplex, [dom])...)*dom.dr1*dom.dz1
        Hij += calc_H0(ϕ[simplex]..., Pi.(simplex, [dom])...)
    end
    return Hij*dom.dr1*dom.dz1
end

"""
    compute_discrete_H(ϕ, dom)

Compute the discrete Heaviside H across the domain, for use in volume integrals. 
"""
function compute_discrete_H(ϕ, dom)
    compute_local_H.(CartesianIndices(ϕ), [ϕ], [dom])
end