export compute_icesurf_δ, compute_icevol_H

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
    outsurf::eltype(ϕ) = 0.0
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
    compute_icesurf_δ(ϕ, dom)

Compute the surface area of the ice, using a discrete Dirac delta function.
Calls [`compute_discrete_δ`](@ref), which is implemented according
to [Min and Gibou 2008](@cite min_robust_2008).
"""
function compute_icesurf_δ(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    SA = 2π*sum(δ .* dom.rgrid)*dom.dr*dom.dz
    return SA
end

"""
    compute_icevol_H(ϕ, dom)

Compute the volume of the ice, using a discrete Heaviside function.
Calls [`compute_discrete_H`](@ref), which is implemented according
to [Min and Gibou 2008](@cite min_robust_2008).
"""
function compute_icevol_H(ϕ, dom)
    H = compute_discrete_H(ϕ, dom)
    vol = 2π*sum(H .* dom.rgrid)*dom.dr*dom.dz
    return vol
end

function compute_icegl_area_weights(ϕ, dom)
    zweights = zeros(eltype(ϕ), dom.nz) 
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
    rweights = zeros(eltype(ϕ), dom.nr) 
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
function compute_icetop_area_weights(ϕ, dom)
    rweights = zeros(eltype(ϕ), dom.nr) 
    for ir in 1:dom.nr-1
        ϕl, ϕr = ϕ[ir:ir+1, dom.nz]
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

Compute the height of ice at each point radially, as well as identify whether it touches system boundary or not.

Returns (`heights`, `bottom_contact`, `top_contact`), where 
- `heights` is a vector of floats
- `bottom_contact`, `top_contact` are vectors of bools
"""
function compute_iceht_bottopcont(ϕ, dom)
    bottom_contact = (ϕ[:,begin] .<= 0)
    top_contact = (ϕ[:,end] .<= 0)

    heights = zeros(eltype(ϕ), dom.nr)
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
                This case is not properly treated." interfaces ir 
            println(ϕ[ir, :] .<= 0)
                
        end
    end
    # if !any(bottom_contact)
    #     display(heat(ϕ, dom))
    # end
    return heights, bottom_contact, top_contact
end


# --------------- Min & Gibou 2008 on surface integrals
# These functions are used only internally by calc_δ0 and calc_H0.

Pij(Pi, Pj, ϕi, ϕj) = Pi*ϕj/(ϕj-ϕi) + Pj*ϕi/(ϕi-ϕj)
Pi(i, j, dom) = [dom.rgrid[i], dom.zgrid[j]]
Pi(ij, dom) = [dom.rgrid[Tuple(ij)[1]], dom.zgrid[Tuple(ij)[2]]]
distance(Pi, Pj) = norm(Pi.-Pj, 2)
area(Pi, Pj, Pk) = abs(Pi[1]*(Pj[2]-Pk[2]) + Pj[1]*(Pk[2]-Pi[2]) + Pk[1]*(Pi[2]-Pj[2]))/2

"""
    calc_δ0(ϕ0, ϕ1, ϕ2, P0, P1, P2)

Implementation of [Min & Gibou 2008, Table 1](@cite min_robust_2008).
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

Implementation of the final expression for a discrete delta in [Min and Gibou 2008](@cite min_robust_2008).
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

Compute the discrete Dirac δ throughout the domain, for use in surface integrals. 
Implements the discrete delta of [Min and Gibou 2008](@cite min_robust_2008).
"""
function compute_discrete_δ(ϕ, dom)
    compute_local_δ.(CartesianIndices(ϕ), [ϕ], [dom])
end

"""
    calc_H0(ϕ0, ϕ1, ϕ2, P0, P1, P2)

Implementation of [Min & Gibou 2008, Table 1](@cite min_robust_2008).
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

Implementation of the final expression for a discrete Heaviside in [Min & Gibou, 2008](@cite min_robust_2008).
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

Implements the discrete Heaviside of [Min & Gibou, 2008](@cite min_robust_2008).
"""
function compute_discrete_H(ϕ, dom)
    compute_local_H.(CartesianIndices(ϕ), [ϕ], [dom])
end