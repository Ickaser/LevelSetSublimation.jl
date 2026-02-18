export identify_Γ, Γ_cells, identify_B 
export reinitialize_ϕ_HCR!, reinitialize_ϕ_HCR
export reinitialize_ϕ_HCR_blindspots!, reinitialize_ϕ_HCR_blindspots

export 𝒢_1st, 𝒢_weno, 𝒢_1st_all, 𝒢_weno_all
export wenodiffs_local
# Functions exported just for the sake of making documentation work


# ---------------- Drawn from Hartmann, 2008, "Constrained reinitialization"

"""
    function identify_Γ(ϕ, dom::Domain)

Identify cells on the sublimation front (interface), returning as a `Matrix::Bool`.
"""
function identify_Γ(ϕ, dom::Domain)
    locs = similar(ϕ, Bool)
    locs .= false
    sg = sign.(ϕ)
    xshift = sg[1:end-1,:] .* sg[2:end,:] # Multiply neighbors in x
    yshift = sg[:,1:end-1] .* sg[:,2:end] # Multiply neighbors in y
    # Any cell with opposite sign of neighbor (or 0) gets added to Γ
    for i in 1:dom.nr-1, j in 1:dom.nz
        if xshift[i,j] <= 0
            locs[i,j] = locs[i+1,j] = true
        end
    end
    for i in 1:dom.nr, j in 1:dom.nz-1
        if yshift[i,j] <= 0
            locs[i,j] = locs[i, j+1] = true
        end
    end
    return locs
end

"""
    Γ_cells(ϕ, dom::Domain) 
    
Compute `findall(identify_Γ(ϕ, dom))`. (That's the whole implementation.)
"""
Γ_cells(ϕ, dom::Domain) = findall(identify_Γ(ϕ, dom))

"""
    identify_B(Γc::Vector{CartesianIndex{2}}, dom::Domain; extra=0)
    identify_B(Γ_field::Matrix{Bool}, dom::Domain)
    identify_B(ϕ::Matrix{Float64}, dom::Domain)
    identify_B(ϕ::Any, dom::Domain)

Return a field of bools identifying the band around the interface.

The width in the band around Γ is specified by the fields `bwr` and `bwz`, 
which represent number of cells in the 𝑟 and 𝑧 directions respectively.
`extra` will tack on extra cells, if in some (but not all) places you need a larger band than the default.
"""
function identify_B(Γc::Vector{CartesianIndex{2}}, dom::Domain; extra=0)
    # nx, ny = size(Γ_field)
    nx = dom.nr
    ny = dom.nz
    B = fill(false, nx, ny)
    # Γc = findall(Γ_field)
    for c in Γc
        ix, iy = Tuple(c)
        xgrab = range(max(1,ix-dom.bwr-extra), min(nx, ix+dom.bwz+extra))
        ygrab = range(max(1,iy-dom.bwr-extra), min(ny, iy+dom.bwz+extra))
        B[xgrab, iy] .= true
        B[ix, ygrab] .= true
    end
    return B
end
function identify_B(Γ_field::Matrix{Bool}, dom::Domain; kwargs...)
    return identify_B(findall(Γ_field), dom; kwargs...)
end
function identify_B(ϕ::Matrix{Float64}, dom::Domain; kwargs...)
    return identify_B(Γ_cells(ϕ, dom), dom;kwargs...)
end
function identify_B(ϕ, dom::Domain;kwargs...)
    return identify_B(Γ_cells(ϕ, dom), dom;kwargs...)
end

"""
    calc_rij_Sij(ϕ, Γ, C, dom::Domain)

Compute rij, neighbors Sij for each cell in `Γ`.

Implementation of eq. 19b from Hartmann 2010, scheme HCR-2
"""
function calc_rij_Sij(ϕ, Γ)
    rij_list = []
    Sij_list = []
    for c in Γ
        pos_neighbors = [CI(1, 0), CI(-1,0), CI(0,1), CI(0,-1)]
        neighbors = [nb for nb in [c].+pos_neighbors if checkbounds(Bool, ϕ, nb)]
        Sij = [nb for nb in neighbors if ϕ[nb]*ϕ[c] <= 0]
        push!(Sij_list, Sij)
        rij = ϕ[c] / sum(ϕ[Sij])
        push!(rij_list, rij)
    end
    return rij_list, Sij_list
end

"""
    reinitialize_ϕ_HCR(ϕ, dom::Domain)

Thin wrapper on `reinitialize_ϕ_HCR!` to avoid mutating.
"""
function reinitialize_ϕ_HCR(ϕ, dom::Domain;kwargs...)
    ϕa = copy(ϕ)
    reinitialize_ϕ_HCR!(ϕa, dom;kwargs...)
    return ϕa
end
"""
    reinitialize_ϕ_HCR(ϕ, dom::Domain)

Thin wrapper on `reinitialize_ϕ_HCR!` to avoid mutating.
"""
function reinitialize_ϕ_HCR_blindspots(ϕ, dom::Domain;kwargs...)
    ϕa = copy(ϕ)
    reinitialize_ϕ_HCR_blindspots!(ϕa, dom;kwargs...)
    return ϕa
end

"""
    sdf_err_L1(ϕ, dom; region=:B)

Compute the L1 norm of the error in the Eikonal equation 
`|∇ϕ| = 1`, on the given region (either `:B` or `:all`).
"""
function sdf_err_L1(ϕ, dom; region=:B)
    if region == :B
        Bf = identify_B(ϕ, dom)
        B = findall(Bf)
        𝒢 = 𝒢_weno.([ϕ], B, [dom])
        return norm(𝒢 .-1, 1) 
    elseif region == :all
        𝒢 = 𝒢_weno_all(ϕ, dom)
        return norm(𝒢 .-1, 1) 
    else
        @error "Bad region to error calc; expect `:B` or `:all`." region
    end
end
"""
    sdf_err_L∞(ϕ, dom; region=:B)

Compute the L∞ norm of the error in the Eikonal equation 
`|∇ϕ| = 1`, on the given region (either `:B` or `:all`).
"""
function sdf_err_L∞(ϕ, dom; region=:B)
    if region == :B
        Bf = identify_B(ϕ, dom)
        B = findall(Bf)
        𝒢 = 𝒢_weno.([ϕ], B, [dom])
    elseif region == :all
        𝒢 = 𝒢_weno_all(ϕ, dom)
    else
        @error "Bad region to error calc; expect `:B` or `:all`." region
    end
    return norm(𝒢 .-1, Inf) 
end

# """
#     calc_err_reg(arr, norm_p=:L∞, region=Colon())

# Compute the norm with given `norm_p` (e.g. `1`, `2`, or `Inf`) of the given array, in provided `region` 
# (which must be a valid set of indices for provided `arr`).

# `region` defaulting to `Colon()` means looking at the full array.
# A `BitArray` or vector of `CartesianIndex`es, as result from e.g. `ϕ .> 0` or `findall(ϕ .> 0)` are valid options.
# """
# function calc_err_reg(arr, norm_p=Inf, region=Colon())
#     # if region == :B
#     #     reg = identify_B(ϕ, dom)
#     # elseif region == :all
#     #     reg = Colon()
#     # else
#     #     @error "Bad region to error calc; expect `:B` or `:all`." region
#     # end

#     # if norm_name == :L∞
#     #     return norm(arr[region], Inf)
#     # elseif norm_name == :L1
#     #     return norm(arr[region], 1)
#     # elseif norm_name == :L2
#     #     return norm(arr[region], 2)
#     # else
#     #     @error "Bad norm to error calc; expect `:B` or `:all`." norm
#     # end
#     norm(arr[region], norm_p)
# end

"""
    reinitialize_ϕ_HCR!(ϕ, dom::Domain; maxsteps = 50, tol=1e-4, err_reg=:B)

Reinitialize `ϕ` throughout the domain.

Implementation of Eq. 22 in Hartmann 2010, scheme HCR-2.

Checks L∞ error (of `|∇ϕ|=1`) against `tol` either in band around interface 
(`err_reg=:B`) or throughout domain (`err_reg=:all`), and ends iteration 
early if tolerance is met.

"""
function reinitialize_ϕ_HCR!(ϕ, dom::Domain; maxsteps = 50, tol=1e-4, err_reg=:B, outside_B = nothing)

    if sum(0 .< extrema(ϕ)) != 1
        @info "Attempted reinit when surface is not in domain. Skipping reinit." extrema(ϕ)
        return nothing
    end
    Γ = Γ_cells(ϕ, dom)
    dx = sqrt(dom.dr*dom.dz) # Geometric mean grid spacing
    Cv = Γ
    F = zeros(size(dom))
    rhs = zeros(size(dom))
    𝒢 = zeros(size(dom))
    S = @. ϕ/sqrt(ϕ^2 + dx^2)


    # Time levels
    dτ = 0.25*dx # Pseudo-time step
    rij_list, Sij_list = calc_rij_Sij(ϕ, Γ)
    signs = sign.(ϕ)

    if err_reg == :B
        region = identify_B(ϕ, dom)
    elseif err_reg == :all
        region = Colon()
    else
        @error "Bad region for error calc; expect `:B` or `:all`." err_reg
    end

    for v in 1:maxsteps
        𝒢 .= 𝒢_weno_all(ϕ, dom; signs=signs)
        if norm(𝒢[region] .- 1, Inf) < tol
            # @info "Early reinit finish. Steps:" v-1
            break
        end
        F .= 0
        rhs .= 0
        for (i,c) in enumerate(Cv)
            # Check for neighbor sign changes, per comment pre Eq. 18
            Sij = Sij_list[i]
            signs_Sij = (ϕ[c] .* ϕ[Sij]) .<= 0
            # If a neighbor no longer has opposite sign, skip this cell
            if sum(signs_Sij) < length(Sij) 
                # @info "Cell skipped" c signs_Sij v
                continue
            end
            # Eq. 21b
            F[c] = (rij_list[i] * sum(ϕ[Sij]) - ϕ[c]) / dx
        end

        # TODO: see Della Rocca and Blanquart 2014 to handle boundaries

        rhs .= S .* (𝒢.-1) .- 0.5F
        if any(isnan.(rhs))
            @warn "NaN in reinit!" findall(isnan.(rhs))
            rhs[isnan.(rhs)] .= 0
        end
        @. ϕ -= rhs * dτ
        # if maximum(abs.(rhs))*dτ < 
        #     @info "Exiting reinit because maxiu"
        # end
    end
    # Outside a band, set to a constant value
    # if err_reg == :B
    #     not_Bf = .~ identify_B(Γ, dom, extra=3)
    #     ϕ[not_Bf] .= sign.(ϕ[not_Bf]) .* (isnothing(outside_B) ? dom.bwfrac*1.5*max(dom.rmax,dom.zmax) : outside_B)
    # end
end

function calc_cosθ_wall(ϕl, ϕnear, dϕdx_m, dϕdx_p, ϕ0sign, dx1) 
    # dw, de, cs = dϕdr_w[cell], dϕdr_e[cell], signs[cell]
    if dϕdx_m > 0 && dϕdx_p > 0 && ϕ0sign > 0
        dϕdη = dϕdx_p
    elseif dϕdx_m > 0 && dϕdx_p > 0 && ϕ0sign < 0
        dϕdη = dϕdx_m
    elseif dϕdx_m < 0 && dϕdx_p < 0 && ϕ0sign > 0
        dϕdη = dϕdx_m
    elseif dϕdx_m < 0 && dϕdx_p < 0 && ϕ0sign < 0
        dϕdη = dϕdx_p
    else
        dϕdη = (abs(dϕdx_m) < abs(dϕdx_p)) ? dϕdx_m : dϕdx_p
    end
    dϕdn = (ϕnear - ϕl)*dx1
    cosθ = dϕdn/norm([dϕdn, dϕdη], 2)
end

"""
    reinitialize_ϕ_HCR_blindspots!(ϕ, dom::Domain; maxsteps = 50, tol=1e-4, err_reg=:B)

Reinitialize `ϕ` throughout the domain.

Implementation of Eq. 22 in Hartmann 2010, scheme HCR-2.
Augmented by (Della Rocca and Blanquart, 2014)'s treatment of contact line blind spots.

Checks L∞ error (of `|∇ϕ|=1`) against `tol` either in band around interface 
(`err_reg=:B`) or throughout domain (`err_reg=:all`), and ends iteration 
early if tolerance is met.

"""
function reinitialize_ϕ_HCR_blindspots!(ϕ, dom::Domain; maxsteps = 50, tol=1e-4, err_reg=:B, outside_B = nothing)

    if sum(0 .< extrema(ϕ)) != 1
        @info "Attempted reinit when surface is not in domain. Skipping reinit." extrema(ϕ)
        return nothing
    end
    Γ = Γ_cells(ϕ, dom)
    dx = sqrt(dom.dr*dom.dz) # Geometric mean grid spacing
    Cv = Γ
    F = zeros(size(dom))
    rhs = zeros(size(dom))
    𝒢 = zeros(size(dom))
    S = @. ϕ/sqrt(ϕ^2 + dx^2)

    # Time levels
    dτ = 0.25*dx # Pseudo-time step
    rij_list, Sij_list = calc_rij_Sij(ϕ, Γ)
    signs = sign.(ϕ)

    dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    # Identify blind walls, per Della Rocca and Blanquart
    # TODO: replace endpoint sign change with internal sign change
    # contactline_south = any(signs[begin:end-1, 1] .* signs[begin+1:end,1] .<= 0)
    # contactline_north = any(signs[begin:end-1, end] .* signs[begin+1:end,end] .<= 0)
    # contactline_west = any(signs[1, begin:end-1] .* signs[1, begin+1:end] .<=0)
    # contactline_east = any(signs[end, begin:end-1] .* signs[end, begin+1:end].<=0)
    blind_south = ((signs[:,1] .* (ϕ[:,2] .- ϕ[:,1])) .>= 0) #.& contactline_south
    blind_north = ((signs[:,end] .* (ϕ[:,end-1] .- ϕ[:,end])) .>= 0) #.& contactline_north
    blind_west =  ((signs[1,:] .* (ϕ[2,:] .- ϕ[1,:])) .>= 0) #.& contactline_west
    blind_east =  ((signs[end,:] .* (ϕ[end-1,:] .- ϕ[end,:])) .>= 0) #.& contactline_east

    blind_south[begin] = blind_south[end] = false
    blind_north[begin] = blind_north[end] = false
    blind_west[begin] = blind_west[end] = false
    blind_east[begin] = blind_east[end] = false

    blind_spots = (sum(blind_south) + sum(blind_north) + sum(blind_west) + sum(blind_east) > 0)
    # if blind_spots
        
    if blind_spots
        @info "Blind spots" blind_south blind_north blind_west blind_east sum(blind_west) sum(blind_east)
        cosθ0_south = map(findall(blind_south)) do i
        c = CI(i,1)
        calc_cosθ_wall(ϕ[c], ϕ[c+CI(0,1)], dϕdr_w[c], dϕdr_e[c], signs[c], dom.dz) 
        end
        cosθ0_north = map(findall(blind_north)) do i
        c = CI(i,dom.nz)
        calc_cosθ_wall(ϕ[c], ϕ[c+CI(0,-1)], dϕdr_w[c], dϕdr_e[c], signs[c], dom.dz) 
        end
        cosθ0_west = map(findall(blind_west)) do i
        c = CI(1,i)
        calc_cosθ_wall(ϕ[c], ϕ[c+CI(1,0)], dϕdz_s[c], dϕdz_n[c], signs[c], dom.dr) 
        end
        cosθ0_east = map(findall(blind_east)) do i
        c = CI(dom.nr,i)
        calc_cosθ_wall(ϕ[c], ϕ[c+CI(-1,0)], dϕdz_s[c], dϕdz_n[c], signs[c], dom.dr) 
        end
    end

    if err_reg == :B
        region = identify_B(ϕ, dom)
    elseif err_reg == :all
        region = Colon()
    else
        @error "Bad region for error calc; expect `:B` or `:all`." err_reg
    end

    for v in 1:maxsteps
        dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n = dϕdx_all = dϕdx_all_WENO(ϕ, dom)
        𝒢 .= 𝒢_weno_all(ϕ, dϕdx_all, dom; signs=signs)
        if norm(𝒢[region] .- 1, Inf) < tol
            # @info "Early reinit finish. Steps:" v-1
            break
        end
        F .= 0
        rhs .= 0
        for (i,c) in enumerate(Cv)
            # Check for neighbor sign changes, per comment pre Eq. 18
            Sij = Sij_list[i]
            signs_Sij = (ϕ[c] .* ϕ[Sij]) .<= 0
            # If a neighbor no longer has opposite sign, skip this cell
            if sum(signs_Sij) < length(Sij) 
                # @info "Cell skipped" c signs_Sij v
                continue
            end
            # Eq. 21b
            F[c] = (rij_list[i] * sum(ϕ[Sij]) - ϕ[c]) / dx
        end


        rhs .= S .* (𝒢.-1) .- 0.5F
        if any(isnan.(rhs))
            @warn "NaN in reinit!" findall(isnan.(rhs))
            rhs[isnan.(rhs)] .= 0
        end
        if blind_spots
            if any(blind_south)
                cosθ_south = map(findall(blind_south)) do i
                c = CI(i,1)
                calc_cosθ_wall(ϕ[c], ϕ[c+CI(0,1)], dϕdr_w[c], dϕdr_e[c], signs[c], dom.dz) 
                end
                @. rhs[blind_south,begin] = 1 * (cosθ0_south - cosθ_south)
            end
            if any(blind_north)
                cosθ_north = map(findall(blind_north)) do i
                c = CI(i,dom.nz)
                calc_cosθ_wall(ϕ[c], ϕ[c+CI(0,-1)], dϕdr_w[c], dϕdr_e[c], signs[c], dom.dz) 
                end
                @. rhs[blind_north,end] = -1 * (cosθ0_north - cosθ_north)
            end
            if any(blind_west)
                cosθ_west = map(findall(blind_west)) do i
                c = CI(1,i)
                calc_cosθ_wall(ϕ[c], ϕ[c+CI(1,0)], dϕdz_s[c], dϕdz_n[c], signs[c], dom.dr) 
                end
                @. rhs[begin,blind_west] = 1 * (cosθ0_west - cosθ_west)
            end
            if any(blind_east)
                cosθ_east = map(findall(blind_east)) do i
                c = CI(dom.nz,i)
                calc_cosθ_wall(ϕ[c], ϕ[c+CI(-1,0)], dϕdz_s[c], dϕdz_n[c], signs[c], dom.dr) 
                end
                @. rhs[end,blind_east] =  -1 * (cosθ0_east - cosθ_east)
            end
        end
        @. ϕ -= rhs * dτ
        # if maximum(abs.(rhs))*dτ < 
        #     @info "Exiting reinit because maxiu"
        # end
    end
    # Outside a band, set to a constant value
    # if err_reg == :B
    #     not_Bf = .~ identify_B(Γ, dom, extra=3)
    #     ϕ[not_Bf] .= sign.(ϕ[not_Bf]) .* (isnothing(outside_B) ? dom.bwfrac*1.5*max(dom.rmax,dom.zmax) : outside_B)
    # end
end

"""
    LD{T}

A "little difference", to make Godunov's scheme in [𝒢_weno](@ref) easier to read.
For `a = LD(x::T)`, 
- `a.p = max(x, 0.0)`
- `a.m = min(x, 0.0)`
"""
struct LD{T} # LD short for Little Difference
    p::T
    m::T
end
LD(x) = LD(max(x, 0.0), min(x, 0.0))
"""
    𝒢_1st(ϕ, i, j, dom::Domain) 
    
Compute the norm of the gradient of `ϕ` at point `i, j` by Godunov's scheme to first-order accuracy.

p. 6830 of Hartmann 2008, 10th page of PDF
"""
function 𝒢_1st(ϕ, i, j, dom::Domain) # p. 6830 of Hartmann, 10th page of PDF
    # pcell = ϕ[i,j]
    if i == 1
        a = LD(0)
        b = LD((ϕ[i+1,j] - ϕ[i,j]) * dom.dr1)
    elseif i == dom.nr
        a = LD((ϕ[i,j] - ϕ[i-1,j]) * dom.dr1)
        b = LD(0)
    else
        a = LD((ϕ[i,j] - ϕ[i-1,j]) * dom.dr1)
        b = LD((ϕ[i+1,j] - ϕ[i,j]) * dom.dr1)
    end
    if j == 1
        c = LD(0)
        d = LD((ϕ[i,j+1] - ϕ[i,j]) * dom.dz1)
    elseif j == dom.nz
        c = LD((ϕ[i,j] - ϕ[i,j-1]) * dom.dz1)
        d = LD(0)
    else
        c = LD((ϕ[i,j] - ϕ[i,j-1]) * dom.dz1)
        d = LD((ϕ[i,j+1] - ϕ[i,j]) * dom.dz1)
    end
    # if ϕ[i,j] >= 0
    #     return sqrt(max(a.p^2, b.m^2) + max(c.p^2, d.m^2))
    # else
    #     return sqrt(max(a.m^2, b.p^2) + max(c.m^2, d.p^2))
    # end
    return 𝒢_loc(a, b, c, d, ϕ[i,j])
end

"""
    𝒢_1st_all(ϕ, dom::Domain)

Compute the norm of the gradient of `ϕ` throughout domain by Godunov's scheme to first-order accuracy.

Internally, calls `𝒢_1st` on all computational cells.
"""
function 𝒢_1st_all(ϕ, dom::Domain)
    return reshape([𝒢_1st(ϕ, i, j, dom) for i in 1:dom.nr, j in 1:dom.nz], dom.nr, dom.nz)
end


"""
    weno_Φ(c, d, e, f)

Return a weighted sum of finite differences for a WENO approximation. 
Defined in §3.2 of [hartmannAccuracyEfficiencyConstrained2009](@cite) .
"""
function weno_Φ(c, d, e, f; ϵ=1e-6)
    s0 = 13*(c-d)^2 + 3*(c-3d)^2
    s1 = 13*(d-e)^2 + 3*(d+e)^2
    s2 = 13*(e-f)^2 + 3*(3*e-f)^2
    α0 = 1/(ϵ+s0)^2
    α1 = 6/(ϵ+s1)^2
    α2 = 3/(ϵ+s2)^2
    ω0 = α0 / (α0 + α1 + α2)
    ω2 = α2 / (α0 + α1 + α2)
    return ω0/3*(c-2*d+e) + (ω2-0.5)/6*(d-2*e + f)
end
"""
    wenodiffs_local(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, dx)

Compute one-sided finite differences, using Jiang and Peng's WENO approximation [jiangWeightedENOSchemes2000](@cite).

A relatively easy-to-read reference is §3.2 of [hartmannAccuracyEfficiencyConstrained2009](@cite) .
"""
function wenodiffs_local(u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, dx)

    dx1 = 1/dx
    dx2 = dx1*dx1

    central_part = ( u_m2 - 8u_m1 + 8u_p1 - u_p2) * dx1 / 12

    dd_m2 = (u_m3 -2u_m2 + u_m1) *dx2
    dd_m1 = (u_m2 -2u_m1 + u_0 ) *dx2
    dd_0  = (u_m1 -2u_0  + u_p1) *dx2
    dd_p1 = (u_0  -2u_p1 + u_p2) *dx2
    dd_p2 = (u_p1 -2u_p2 + u_p3) *dx2


    du_l = central_part - dx*weno_Φ(dd_m2, dd_m1, dd_0, dd_p1)
    du_r = central_part + dx*weno_Φ(dd_p2, dd_p1, dd_0, dd_m1)
    
    return du_l, du_r
end

function wenodiffs_row(u, dx)
    dx112 = 1/dx/12
    dx2 = dx^-2

    n = size(u, 1)

    u_ex = similar(u, n + 6) # Extrapolate three points beyond on each side
    u_ex[begin+3:end-3] .= u

    # u_ex[3:-1:1] = extrap_ϕ_mat_quad * u[begin:begin+2] # Quadratic extrapolation, left
    # u_ex[end-2:end] = extrap_ϕ_mat_quad * u[end:-1:end-2] # Quadratic extrapolation, right

    # u_ex[3:-1:1] = extrap_ϕ_mat_lin * u[begin:begin+1] # Linear extrapolation, left
    # u_ex[end-2:end] = extrap_ϕ_mat_lin * u[end:-1:end-1] # Linear extrapolation, right

    u_ex[3:-1:1] .= u[begin] # Constant extrapolation, left
    u_ex[end-2:end] .= u[end] # Constant extrapolation, right

    i0 = (1:n) .+ 3
    central_part = @. (u_ex[i0.-2] - 8u_ex[i0.-1] + 8u_ex[i0.+1] - u_ex[i0.+2]) * dx112

    cdiffs = dx2 * (u_ex[begin:end-2] .- 2u_ex[begin+1:end-1] .+ u_ex[begin+2:end])

    i0 = (1:n) .+ 2
    # left_Φ = weno_Φ.(cdiffs[i0.-2], cdiffs[i0.-1], cdiffs[i0], cdiffs[i0.+1])
    # right_Φ = weno_Φ.(cdiffs[i0.+2], cdiffs[i0.+1], cdiffs[i0], cdiffs[i0.-1])
    # du_l = central_part .- dx*left_Φ
    # du_r = central_part .+ dx*right_Φ
    du_l = central_part .- dx*weno_Φ.(cdiffs[i0.-2], cdiffs[i0.-1], cdiffs[i0], cdiffs[i0.+1])
    du_r = central_part .+ dx*weno_Φ.(cdiffs[i0.+2], cdiffs[i0.+1], cdiffs[i0], cdiffs[i0.-1])

    return du_l, du_r
end

"""
    dϕdx_all_WENO(ϕ, dom)

Compute (`∂ϕ/∂r` west, `∂ϕ/∂r` east, `∂ϕ/∂z` south, `∂ϕ/∂z` north) using WENO derivatives.

Implemented by computing WENO derivatives for each cell separately, which is a little wasteful.
Beyond the boundaries of domain, ϕ is extrapolated according to `get_or_extrapolate_ϕ`.
"""
function dϕdx_all_WENO(ϕ, dom)
    # dϕdr_w = zeros(size(dom))
    # dϕdr_e = zeros(size(dom))
    # dϕdz_s = zeros(size(dom))
    # dϕdz_n = zeros(size(dom))
    dϕdr_w = similar(ϕ)
    dϕdr_e = similar(ϕ)
    dϕdz_s = similar(ϕ)
    dϕdz_n = similar(ϕ)
    # dϕdr_w .= 0
    # dϕdr_e .= 0
    # dϕdz_s .= 0
    # dϕdz_n .= 0

    for i in axes(dϕdr_e, 2)
        dϕdr_w[:,i], dϕdr_e[:,i] = wenodiffs_row(ϕ[:,i], dom.dr)
    end
    for i in axes(dϕdz_n, 1)
        dϕdz_s[i,:], dϕdz_n[i,:] = wenodiffs_row(ϕ[i,:], dom.dz)
    end
    return dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n
end

"""
    𝒢_weno(ϕ, ir::Int, iz::Int, dom::Domain)
    𝒢_weno(ϕ, ind::CartesianIndex{2}, dom::Domain)
    𝒢_weno(ϕ, ir::Int, iz::Int, ϕ0sign, dom::Domain)
    𝒢_weno(ϕ, ind::CartesianIndex{2}, ϕ0sign, dom::Domain)

Compute the norm of the gradient by Godunov's scheme with WENO differences (`wenodiffs_local`).

If supplied, `ϕ0sign` is used in Godunov's scheme, rather than the current sign of ϕ.
Described in [hartmannAccuracyEfficiencyConstrained2009](@cite), eq. 6 to eq. 9.
Boundary treatment of ghost cells handled by `get_or_extrapolate_ϕ`.
"""
function 𝒢_weno(ϕ, ir::Int, iz::Int, ϕ0sign, dom::Domain)
    return 𝒢_weno(ϕ, CI(ir, iz), ϕ0sign, dom)
end
function 𝒢_weno(ϕ, ir::Int, iz::Int, dom::Domain)
    return 𝒢_weno(ϕ, CI(ir, iz), sign(ϕ[ir,iz]), dom)
end
function 𝒢_weno(ϕ, ind::CartesianIndex{2}, dom::Domain)
    return 𝒢_weno(ϕ, ind, sign(ϕ[ind]), dom)
end

function 𝒢_weno(ϕ, ind::CartesianIndex{2}, ϕ0sign, dom::Domain)
    rstenc = [CI(ir, 0) for ir in -3:3]
    zstenc = [CI( 0,iz) for iz in -3:3]

    ar_, br_ = wenodiffs_local(get_or_extrapolate_ϕ(ϕ, ind, rstenc)..., dom.dr)
    az_, bz_ = wenodiffs_local(get_or_extrapolate_ϕ(ϕ, ind, zstenc)..., dom.dz)

    return 𝒢_loc(ar_, br_, az_, bz_, ϕ0sign)
end

function 𝒢_loc(ar, br, az, bz, ϕloc)
    if ϕloc >= 0
        return sqrt(max(ar, -br, 0)^2 + max(az, -bz, 0)^2)
    else
        return sqrt(max(-ar, br, 0)^2 + max(-az, bz, 0)^2)
    end
end
    

"""
    𝒢_weno_all_old(ϕ, dom::Domain)

Compute the norm of the gradient of `ϕ` throughout domain by Godunov's scheme with WENO derivatives.

Internally, calls `𝒢_weno` on all computational cells.
"""
function 𝒢_weno_all_old(ϕ, dom::Domain)
    # return reshape([𝒢_weno(ϕ, i, j, dom) for i in 1:dom.nr, j in 1:dom.nz], dom.nr, dom.nz)
    return 𝒢_weno.([ϕ], CartesianIndices(ϕ), [dom])
end

"""
    𝒢_weno_all(ϕ, [dϕdx_all,] dom::Domain; signs=nothing)

Compute the norm of the gradient of `ϕ` throughout domain by Godunov's scheme with WENO derivatives.
If signs of ϕ should refer to some prior state, can be provided.
If dϕdx_all is provided, then is not recomputed internally.

Internally, calls `𝒢_weno` on all computational cells.
"""
function 𝒢_weno_all(ϕ, dom::Domain; signs=nothing)
    dϕdx_all = dϕdx_all_WENO(ϕ, dom)
    return 𝒢_weno_all(ϕ, dϕdx_all, dom; signs=signs)
end
function 𝒢_weno_all(ϕ, dϕdx_all, dom::Domain; signs=nothing)
    if isnothing(signs)
        𝒢 = 𝒢_loc.(dϕdx_all..., ϕ)
    else
        𝒢 = 𝒢_loc.(dϕdx_all..., signs)
    end
    return 𝒢
end

"""
    dϕdx_all_WENO_loc(ϕ, dom)

Compute (`∂ϕ/∂r` west, `∂ϕ/∂r` east, `∂ϕ/∂z` south, `∂ϕ/∂z` north) using WENO derivatives.

Implemented by computing WENO derivatives for each cell separately, which is a little wasteful.
Beyond the boundaries of domain, ϕ is extrapolated according to `get_or_extrapolate_ϕ`.
"""
function dϕdx_all_WENO_loc(ϕ, dom)
    # --- Compute ϕ derivatives with WENO
    # indmin = CI(1, 1)
    # indmax = CI(dom.nr, dom.nz)
    rstenc = [CI(i, 0) for i in -3:3]
    zstenc = [CI(0, i) for i in -3:3]

    dϕdr_w = zeros(Float64, size(dom))
    dϕdr_e = zeros(Float64, size(dom))
    dϕdz_s = zeros(Float64, size(dom))
    dϕdz_n = zeros(Float64, size(dom))
    for ind in CartesianIndices(ϕ)
        # rst = max.(min.([ind].+rshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
        # zst = max.(min.([ind].+zshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
        # dϕdr_w[ind], dϕdr_e[ind] = wenodiffs_local(ϕ[rst]..., dom.dr)
        # dϕdz_s[ind], dϕdz_n[ind] = wenodiffs_local(ϕ[zst]..., dom.dz)

        dϕdr_w[ind], dϕdr_e[ind] = wenodiffs_local(get_or_extrapolate_ϕ(ϕ, ind, rstenc)..., dom.dr)
        dϕdz_s[ind], dϕdz_n[ind] = wenodiffs_local(get_or_extrapolate_ϕ(ϕ, ind, zstenc)..., dom.dz)
    end
    return dϕdr_w, dϕdr_e, dϕdz_s, dϕdz_n
end

const extrap_ϕ_mat_quad   = [3 -3 1 ; 6 -8 3 ; 10 -15 6] # quadratic
const extrap_ϕ_mat_lin = [2 -1; 3 -2; 4 -3] # linear
""" 
    get_or_extrapolate_ϕ(ϕ, ind, stencil)

Either retrieve `ϕ` at `ind` on the given `stencil`, or extrapolate and return domain + extrapolated values.

This extrapolation is on a uniform grid, using a quadratic extrapolant.
Define Lagrange interpolant with last three points in domain, then 
extrapolate outside the domain to make ghost points.
Ghost points are then a linear combination of points from domain.
A matrix form just keeps this compact and efficient;
`extrap_ϕ_mat` in the source is just the matrix representing this extrapolation.
""" 
function get_or_extrapolate_ϕ(ϕ, ind, stencil)
        # See how much of the stencil is inside the domain
        inside = checkbounds.(Bool, [ϕ], [ind].+stencil)
        # If all of the stencil is in the domain, no extrapolation needed
        if findlast(inside) - findfirst(inside) == length(stencil) - 1
            # @info "interior" [ind] .+ stencil
            return ϕ[[ind].+stencil]
        end
        # If some of the stencil is outside domain,
        ϕis = similar(ϕ, size(stencil))
        i1 = findfirst(inside)
        il = findlast(inside)
        # use as much of the domain as possible,
        ϕis[i1:il] .= ϕ[[ind].+stencil[i1:il]]

        # then extrapolate for as many points as necessary
        if firstindex(ϕis) != i1
            # ϕis[i1-1:-1:begin] = (extrap_ϕ_mat_quad[1:i1-1,:] * ϕis[i1:i1+2]) # Quadratic extrapolation
            # ϕis[i1-1:-1:begin] = (extrap_ϕ_mat_lin[1:i1-1,:] * ϕis[i1:i1+1]) # LInear extrapolation
            ϕis[i1-1:-1:begin] .= ϕis[i1] # Constant extrapolation
        elseif lastindex(ϕis) != il
            # ϕis[il+1:end] = extrap_ϕ_mat_quad[1:lastindex(ϕis)-il,:] * ϕis[il:-1:il-2] # Quadratic extrapolation
            # ϕis[il+1:end] = extrap_ϕ_mat_lin[1:lastindex(ϕis)-il,:] * ϕis[il:-1:il-1] # LInear extrapolation
            ϕis[il+1:end] .= ϕis[il] # Constant extrapolation
        end

        return ϕis
end
