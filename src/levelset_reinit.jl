export identify_Γ, Γ_cells, identify_B 
export reinitialize_ϕ_HCR!, reinitialize_ϕ_HCR

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
    identify_B(Γc::Vector{CartesianIndex{2}}, dom::Domain)
    identify_B(Γ_field::Matrix{Bool}, dom::Domain)
    identify_B(ϕ::Matrix{Float64}, dom::Domain)
    identify_B(ϕ::Any, dom::Domain)

Return a field of bools identifying the band around the interface.

The width in the band around Γ is specified by the fields `bwr` and `bwz`, 
which represent number of cells in the 𝑟 and 𝑧 directions respectively.
"""
function identify_B(Γc::Vector{CartesianIndex{2}}, dom::Domain)
    # nx, ny = size(Γ_field)
    nx = dom.nr
    ny = dom.nz
    B = fill(false, nx, ny)
    # Γc = findall(Γ_field)
    for c in Γc
        ix, iy = Tuple(c)
        xgrab = range(max(1,ix-dom.bwr), min(nx, ix+dom.bwz))
        ygrab = range(max(1,iy-dom.bwr), min(ny, iy+dom.bwz))
        B[xgrab, iy] .= true
        B[ix, ygrab] .= true
    end
    return B
end
function identify_B(Γ_field::Matrix{Bool}, dom::Domain)
    return identify_B(findall(Γ_field), dom)
end
function identify_B(ϕ::Matrix{Float64}, dom::Domain)
    return identify_B(Γ_cells(ϕ, dom), dom)
end
function identify_B(ϕ, dom::Domain)
    return identify_B(Γ_cells(ϕ, dom), dom)
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
function reinitialize_ϕ_HCR(ϕ, dom::Domain)
    ϕa = copy(ϕ)
    reinitialize_ϕ_HCR!(ϕa, dom)
    return ϕa
end

function sdf_err_L1(ϕ, dom; region=:B)
    if region == :B
        Bf = identify_B(ϕ, dom)
        B = findall(Bf)
        𝒢 = 𝒢_weno.([ϕ], B, [dom])
    elseif region == :all
        𝒢 = 𝒢_weno_all(ϕ, dom)
    else
        @error "Bad region to error calc; expect `:B` or `:all`." region
    end
    err = sum(abs.(𝒢 .-1)) / length(B)
end
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
    err = maximum(abs.(𝒢 .-1)) 
end

"""
    reinitialize_ϕ_HCR2!(ϕ, dom::Domain; maxsteps = 20, tol=1e-4, err_reg=:B)

Reinitialize `ϕ` throughout the domain.

Implementation of Eq. 22 in Hartmann 2010, scheme HCR-2.

Checks L∞ error against `tol` either in band around interface (`err_reg=:B`) or throughout domain (`err_reg=:all`), and ends iteration early

TODO: switch to Eq. 23 to minimize allocations? Can eliminate F, rhs that way

"""
function reinitialize_ϕ_HCR!(ϕ, dom::Domain; maxsteps = 50, tol=1e-4, err_reg=:B)

    if sum(0 .< extrema(ϕ)) != 1
        @info "Attempted reinit when surface is not in domain. Skipping reinit." extrema(ϕ)
        return nothing
    end
    Γ = Γ_cells(ϕ, dom)
    dx = sqrt(dom.dr*dom.dz) # Geometric mean grid spacing
    Cv = Γ
    F = zeros(size(dom))
    rhs = zeros(size(dom))
    S = @. ϕ/sqrt(ϕ^2 + dx^2)
    # Time levels
    dτ = 0.25*dx # Pseudo-time step
    rij_list, Sij_list = calc_rij_Sij(ϕ, Γ)
    # sdf_err_L1 = 
    for v in 1:maxsteps
        # if sdf_err_L1(ϕ, dom) < tol
        if sdf_err_L∞(ϕ, dom, region=err_reg) < tol
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
        rhs .= S .* (𝒢_weno.([ϕ], CartesianIndices(ϕ), [dom]) .- 1) .- 0.5F
        ϕ .-= rhs .* dτ
    end
end


"""
    LD{T}

    A "little difference", to make Godunov's scheme in [𝒢_weno](@ref) easier to read.
For a = LD(x::T), 
- a.p = max(x, 0)
- a.m = min(x, 0)
"""
struct LD{T} # LD short for Little Difference
    p::T
    m::T
end
LD(x) = LD(max(x, 0), min(x, 0))
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
    if ϕ[i,j] >= 0
        return sqrt(max(a.p^2, b.m^2) + max(c.p^2, d.m^2))
    else
        return sqrt(max(a.m^2, b.p^2) + max(c.m^2, d.p^2))
    end
end

"""
    𝒢_1st_all(ϕ, dom::Domain)

Compute the norm of the gradient of `ϕ` throughout domain by Godunov's scheme to first-order accuracy.

Internally, calls [𝒢_1st](@ref) on all computational cells.
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
"""
    𝒢_weno(ϕ, ir::Int, iz::Int, dom::Domain)
    𝒢_weno(ϕ, ind::CartesianIndex{2}, dom::Domain)

Compute the norm of the gradient by Godunov's scheme with WENO differences ([`wenodiffs_local`](@ref)).

Described in [hartmannAccuracyEfficiencyConstrained2009](@cite), eq. 6 to eq. 9.
Let all ghost cells equal the function value at boundary; I think this is equivalent to using homogeneous Neumann boundaries.
"""
function 𝒢_weno(ϕ, ir::Int, iz::Int, dom::Domain)
    return 𝒢_weno(ϕ, CI(ir, iz), dom)
end


function 𝒢_weno(ϕ, ind::CartesianIndex{2}, dom::Domain)
    indmin = CI(1, 1)
    indmax = CI(dom.nr, dom.nz)
    rshift = [CI(ir, 0) for ir in -3:3]
    zshift = [CI( 0,iz) for iz in -3:3]
    rst = max.([indmin], min.([indmax], [ind].+rshift)) # If stencil falls partly outside domain,  
    zst = max.([indmin], min.([indmax], [ind].+zshift)) # repeat the boundary cell

    ar_, br_ = wenodiffs_local(ϕ[rst]..., dom.dr)
    az_, bz_ = wenodiffs_local(ϕ[zst]..., dom.dz)

    ar = LD(ar_)
    br = LD(br_)
    az = LD(az_)
    bz = LD(bz_)

    if ϕ[ind] >= 0
        return sqrt(max(ar.p^2, br.m^2) + max(az.p^2, bz.m^2))
    else
        return sqrt(max(ar.m^2, br.p^2) + max(az.m^2, bz.p^2))
    end
end

"""
    𝒢_weno_all(ϕ, dom::Domain)

Compute the norm of the gradient of `ϕ` throughout domain by Godunov's scheme with WENO derivatives.

Internally, calls `𝒢_weno` on all computational cells.
"""
function 𝒢_weno_all(ϕ, dom::Domain)
    return reshape([𝒢_weno(ϕ, i, j, dom) for i in 1:dom.nr, j in 1:dom.nz], dom.nr, dom.nz)
    # return 𝒢_weno([ϕ], CartesianIndices(ϕ), [dom])
end