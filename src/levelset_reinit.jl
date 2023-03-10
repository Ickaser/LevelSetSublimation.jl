export identify_Γ, Γ_cells, identify_B, plot_RC, 𝒢_all
export reinitialize_ϕ, reinitialize_ϕ!, reinitialize_ϕ_all!
export reinitialize_ϕ_HCR!, reinitialize_ϕ_HCR

# Functions exported just for the sake of making documentation work
export update_ϕ_in_Γ!
export calc_dϕdr_sdf, calc_dϕdz_sdf, identify_regions_RC
export 𝒢_1st, 𝒢_weno, 𝒢_1st_all, 𝒢_weno_all


# ---------------- Drawn from Hartmann, 2008, "Constrained reinitialization"

"""
    function identify_Γ(ϕ, dom::Domain)

Identify cells on the sublimation front (interface), returning as a `Matrix::Bool`.
"""
function identify_Γ(ϕ, dom::Domain)
    locs = similar(ϕ, Bool)
    locs .= false
    sg = sign.(ϕ)
    xshift = sg[1:end-1,:] .* sg[2:end,:]
    yshift = sg[:,1:end-1] .* sg[:,2:end]
    
    for i in 1:dom.nr-1, j in 1:dom.nz
        # if locs[i,j]
        #     continue
        # end
        if xshift[i,j] <= 0
            locs[i,j] = locs[i+1,j] = true
        end
    end
    for i in 1:dom.nr, j in 1:dom.nz-1
        # if locs[i,j]
        #     continue
        # end
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
    function calc_curvature(ϕ, dom::Domain)

Not used explicitly at present, but useful for debugging.
"""
function calc_curvature(ϕ, dom::Domain)
    dx2 = dom.dr2
    dy2 = dom.dz2
    # dx2 = 1/dx^2
    # dy2 = 1/dy^2
    nx = dom.nr
    ny = dom.nz
    # nx, ny = size(ϕ)
    ∇2ϕ = similar(ϕ)
    
    # Second order everywhere (upwind at edges)

    for ix in 2:nx-1
        ∇2ϕ[ix,:] = @. (ϕ[ix+1,:] - 2ϕ[ix,:] + ϕ[ix-1,:] )*dx2
    end
    ∇2ϕ[1,:] = ∇2ϕ[2,:] # Uses the same stencil, unfortunately
    ∇2ϕ[end,:] = ∇2ϕ[end-1,:] # Uses the same stencil, unfortunately

    # y portion, need to be slightly more careful about reusing stencils
    ∇2ϕ[:,1] += @. (ϕ[:,3] - 2ϕ[:,2] + ϕ[:,1])*dy2
    ∇2ϕ[:,end] += @. (ϕ[:,end] - 2ϕ[:,end-1] + ϕ[:,end-2])*dy2
    for iy in 2:ny-1
        ∇2ϕ[:,iy] += @. (ϕ[:,iy+1] -2ϕ[:,iy] + ϕ[:,iy-1] )*dy2
    end
            
    return -∇2ϕ
end

"""
    function identify_regions_RC(ϕ, Γ, dom::Domain)

Takes full level set field ϕ, list of front cells Γ, and domain.
Computes curvature (or at least something proportional to it) at all locations Γ, then compares against sign of ϕ to assign to R or C
"""
function identify_regions_RC(ϕ, Γ, dom::Domain)
    # dx2 = 1/dx^2
    # dy2 = 1/dy^2
    dx2 = dom.dr2
    dy2 = dom.dz2
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    numcells = length(Γ)
    CC = fill(0.0, numcells) # Note: Curvature = -∇^2(ϕ)
    R = Vector{CartesianIndex{2}}()
    C = Vector{CartesianIndex{2}}()
    for ic in 1:numcells
        cell = Γ[ic]
        ix, iy = Tuple(cell)
        if ix == 1
            CC[ic] -= (ϕ[3,iy] - 2ϕ[2,iy] + ϕ[1,iy])*dx2
        elseif ix == nx
            CC[ic] -= (ϕ[nx,iy] - 2ϕ[nx-1,iy] + ϕ[nx-2,iy])*dx2
        else
            CC[ic] -= (ϕ[ix+1,iy] - 2ϕ[ix,iy] + ϕ[ix-1,iy])*dx2
        end
        if iy == 1
            CC[ic] -= (ϕ[ix,3] - 2ϕ[ix,2] + ϕ[ix,1])*dy2
        elseif iy == ny
            CC[ic] -= (ϕ[ix,ny] - 2ϕ[ix,ny-1] + ϕ[ix,ny-2])*dy2
        else
            CC[ic] -= (ϕ[ix,iy+1] - 2ϕ[ix,iy] + ϕ[ix,iy-1])*dx2
        end

        CC = round.(CC, digits=7)

        if CC[ic]*ϕ[cell] < 0 || (ϕ[cell] < 0 && CC[ic] == 0 )
            push!(C, cell)
            # println("C: c=$((ix, iy)), CC[c] = $(CC[ic]), ϕ[c] = $(ϕ[cell])")
        else
            push!(R, cell)
            # println("R: c=$((ix, iy)), CC[c] = $(CC[ic]), ϕ[c] = $(ϕ[cell])")
        end
    end
    return R, C

end

"""
    function plot_RC(ϕ, dom::Domain)

Mark cells in ℛ with black, cells in 𝒞 with white, in mutating fashion.
"""
function plot_RC(ϕ, dom::Domain)
    R, C = identify_regions_RC(ϕ, Γ_cells(ϕ, dom), dom)
    Rr = [rgrid[Tuple(c)[1]] for c in R]
    Rz = [zgrid[Tuple(c)[2]] for c in R]
    Cr = [rgrid[Tuple(c)[1]] for c in C]
    Cz = [zgrid[Tuple(c)[2]] for c in C]
    scatter!(Rr, Rz, c=:black)
    scatter!(Cr, Cz, c=:white)
    # plot_RC(RC, nr, nz)
end

"""
    identify_B(Γc::Vector{CartesianIndex{2}}, dom::Domain)
    identify_B(Γ_field::Matrix{Bool}, dom::Domain)
    identify_B(ϕ::Matrix{Float64}, dom::Domain)

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

"""
    calc_dϕdr_sdf(ϕ, Γf, i, j, dom::Domain)

Take a derivative in 𝑟 inside Γ, for computing signed distance function.

This is equivalent to part of Eq. 21 in Hartmann 2008.
"""
function calc_dϕdr_sdf(ϕ, Γf, i, j, dom::Domain)
    # Fancy conditions for near coalescence: ignored for now
    # TODO: fill these out for real. If A && B is true, jp = j, likewise for jm = j or something
    # A = true
    # B = true
    if i == dom.nr
        ip = i
        im = (Γf[i-1,j] ? i-1 : i)
    elseif i == 1
        ip = (Γf[i+1,j] ? i+1 : i)
        im = i
    else
        ip = (Γf[i+1,j] ? i+1 : i)
        im = (Γf[i-1,j] ? i-1 : i)
    end
    num = ϕ[ip,j] - ϕ[im,j]
    den = max(dom.rgrid[ip] - dom.rgrid[im], .001*dom.dr)
    return num/den
end

"""
    calc_dϕdr_sdf(ϕ, Γf, i, j, dom::Domain)

Take a derivative in 𝑧 inside Γ, for computing signed distance function.

This is equivalent to part of Eq. 21 in Hartmann 2008.
"""
function calc_dϕdz_sdf(ϕ, Γf, i, j, dom::Domain)
    # Fancy conditions for near coalescence: ignored for now
    # TODO: fill these out for real. If A && B is true, jp = j, likewise for jm = j or something
    # A = true
    # B = true
    if j == dom.nz
        jp = j
        jm = (Γf[i,j-1] ? j-1 : j)
    elseif j == 1
        jp = (Γf[i,j+1] ? j+1 : j)
        jm = j
    else
        jp = (Γf[i,j+1] ? j+1 : j)
        jm = (Γf[i,j-1] ? j-1 : j)
    end
    num = ϕ[i,jp] - ϕ[i,jm]
    den = max(dom.zgrid[jp] - dom.zgrid[jm], .001*dom.dz)
    return num/den
end

function calc_dij_R!(d, ϕ, Γf, R, dom::Domain)
    # nr, nz = size(ϕ)
    for c in R
        if ϕ[c]==0
            d[c] = 0
            continue
        end
        i, j = Tuple(c)
        den = hypot(calc_dϕdr_sdf(ϕ, Γf, i, j, dom), calc_dϕdz_sdf(ϕ, Γf, i, j, dom))
        d[c] = ϕ[c] / den
    end
end
function calc_dij_C!(d, ϕ, C, dom::Domain)
    # nr, nz = size(ϕ)
    for c in C
        if(ϕ[c]==0)
            d[c] = 0
            continue
        end
        pos_neighbors = [CI(1, 0), CI(-1,0), CI(0,1), CI(0,-1)]
        neighbors = [nb for nb in [c].+pos_neighbors if checkbounds(Bool, ϕ, nb)]

        if length(Sij) > 0
            num = sum([d[nb] for nb in Sij])
            den = sum([ϕ[nb] for nb in Sij])
            d[c] =  ϕ[c] * num / den
        else
            @warn "Length of Sij is $(length(Sij))! Not updating value."
            d[c] = ϕ[c]
            continue
        end
    end
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

function sdf_err_L1(ϕ, dom)
    Bf = identify_B(ϕ, dom)
    B = findall(Bf)
    𝒢 = 𝒢_weno.([ϕ], B, [dom])
    err = sum(abs.(𝒢 .-1)) / length(B)
end
function sdf_err_L∞(ϕ, dom)
    Bf = identify_B(ϕ, dom)
    B = findall(Bf)
    𝒢 = 𝒢_weno.([ϕ], B, [dom])
    err = maximum(abs.(𝒢 .-1)) 
end

"""
    reinitialize_ϕ_HCR2!(ϕ, dom::Domain; maxsteps = 20)

Reinitialize `ϕ` throughout the domain.

Implementation of Eq. 22 in Hartmann 2010, scheme HCR-2.

TODO: switch to Eq. 23 to minimize allocations? Can eliminate F, rhs that way

"""
function reinitialize_ϕ_HCR!(ϕ, dom::Domain; maxsteps = 20, tol=1e-4)
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
        if sdf_err_L1(ϕ, dom) < tol
            @info "End reinit early" sdf_err_L1(ϕ, dom) v
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
                continue
            end
            # Eq. 21b
            F[c] = (rij_list[i] * sum(ϕ[Sij]) - ϕ[c]) / dx
            # @info "F" c F[c] rij_list[i]*sum(ϕ[Sij])
        end
        # for c in CartesianIndices(ϕ)
        #     𝒢 = 𝒢_weno(ϕ, c, dom)
        #     rhs[c] = dτ * (S[c]*(𝒢 - 1) - 0.5F[c])
        # end
        # 𝒢 = 𝒢_weno.([ϕ], CartesianIndices(ϕ), [dom]) .- 1
        # rhs .= S .* 𝒢 .- 0.5F
        rhs .= S .* (𝒢_weno.([ϕ], CartesianIndices(ϕ), [dom]) .- 1) .- 0.5F
        # @info "step" S F 𝒢 rhs 
        # @info "Timestep" v F
        ϕ .-= rhs .* dτ
    end
end


"""
    update_ϕ_in_Γ!(ϕ, dom::Domain)

Reinitialize the interface cells of `ϕ`.

This is the scheme CR-2 in Hartmann 2008 (note the published erratum to that article, which amends ℛ and 𝒞).
"""
function update_ϕ_in_Γ!(ϕ, dom::Domain)
    Γf = identify_Γ(ϕ, dom)
    Γc = findall(Γf)
    RC = identify_regions_RC(ϕ, Γc, dom)
    dl = fill(0.0, dom.nr, dom.nz)
    calc_dij_R!(dl, ϕ, Γf, Γc, dom)
    calc_dij_C!(dl, ϕ, RC[2], dom)
    for c in Γc
        ϕ[c] = dl[c]
    end
end

"""
    LD{T}

    A "little difference", to make Godunov's scheme in [𝒢](@ref) easier to read.
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
    reinitialize_ϕ!(ϕ, dom::Domain, tf=1.0; alg=BS3(), outside_B = 1)

Reinitialize the signed distance function `ϕ`.

Carried out in place.  

This relies on `update_ϕ_in_Γ`, which implements CR-2 from Hartmann 2008,
then in a band ℬ around the interface, solves a reinitialization PDE
using a first-order or WENO spatial scheme with time integration given by `alg`.
"""
function reinitialize_ϕ!(ϕ, dom::Domain, tf=100.0; alg=BS3())

    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Bf = identify_B(Γ, dom)
    BnΓ = findall(Bf .⊻ Γf)
    ΩnB = findall(fill(true, dom.nr, dom.nz) .⊻ Bf)

    outside_B = 1.5*dom.bwfrac*max(dom.rmax, dom.zmax)
    update_ϕ_in_Γ!(ϕ, dom)

    sarr = sign.(ϕ)
    Γ = Γ_cells(ϕ, dom)

    
    ϕ_ode = ϕ[BnΓ]
    cached = copy(ϕ)
    function sub_rhs(du, u, p, t) 
        cached[BnΓ] .= u
        for (i, c) in enumerate(BnΓ)
            # du[i] = sarr[c] * (1-𝒢_1st(cached, Tuple(c)..., dom))
            du[i] = sarr[c] * (1-𝒢_weno(cached, Tuple(c)..., dom))
        end
        return du
    end
    tspan = (0.0, tf)
    prob = ODEProblem(sub_rhs, ϕ_ode, tspan)
    sol = solve(prob, alg, dt = 1.0; callback=TerminateSteadyState(1e-4, 1e-4))
    ϕ[BnΓ] .= sol[end]

    # @info "Reinitialization time" sol.t[end]

    ϕ[ΩnB] .= sarr[ΩnB] .* outside_B

    nothing
end

"""
    reinitialize_ϕ(ϕ, dom::Domain, tf=1.0; alg=BS3(), outside_B = 1)

Reinitialize the signed distance function `ϕ`, returning a new array.

Simply makes a copy of `ϕ` and calls `reinitialize_ϕ!`.
"""
function reinitialize_ϕ(ϕ, dom::Domain, tf=1.0; alg = BS3())
    ϕ1 = copy(ϕ)
    reinitialize_ϕ!(ϕ1, dom, tf; alg=alg)
    ϕ1
end


"""
    reinitialize_ϕ_all!(ϕ, dom::Domain, tf=1.0; alg=BS3(), outside_B = 1)

Reinitialize the signed distance function `ϕ`.

Carried out in place.  

This relies on `update_ϕ_in_Γ`, which implements CR-2 from Hartmann 2008,
then everywhere else, solves a reinitialization PDE
using a WENO spatial scheme with time integration given by `alg`.
"""
function reinitialize_ϕ_all!(ϕ, dom::Domain, tf=100.0; alg=BS3())
    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    sarr = sign.(ϕ)

    update_ϕ_in_Γ!(ϕ, dom)

    function sub_rhs(du, u, p, t) 
        ϕl = reshape(u, dom.nr, dom.nz)
        dϕ = sarr .* (1 .- 𝒢_weno_all(ϕl, dom))
        dϕ[Γ] .= 0.0
        du .= reshape(dϕ, :)
        return du
    end
    tspan = (0.0, tf)
    ϕ_flat = reshape(ϕ, :)
    prob = ODEProblem(sub_rhs, ϕ_flat, tspan)
    sol = solve(prob, alg, dt = 1.0; callback=TerminateSteadyState(1e-4, 1e-4))
    ϕ .= reshape(sol[end], dom.nr, dom.nz)

    # @info "Reinitialization time" sol.t[end]

    nothing
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

Compute the norm of the gradient by Godunov's scheme with WENO differences ([wenodiffs_local](@ref)).

Described in [hartmannAccuracyEfficiencyConstrained2009](@cite), eq. 6 to eq. 9.
Let all ghost cells equal the function value at boundary; I think this is equivalent to using homogeneous Neumann boundaries.
"""
function 𝒢_weno(ϕ, ir::Int, iz::Int, dom::Domain)
    irs = max.(1, min.(dom.nr, ir-3:ir+3)) # Pad with boundary values
    izs = max.(1, min.(dom.nz, iz-3:iz+3))

    ar, br = wenodiffs_local(ϕ[irs, iz]..., dom.dr)
    az, bz = wenodiffs_local(ϕ[ir, izs]..., dom.dz)

    ar = LD(ar)
    br = LD(br)
    az = LD(az)
    bz = LD(bz)

    if ϕ[ir,iz] >= 0
        return sqrt(max(ar.p^2, br.m^2) + max(az.p^2, bz.m^2))
    else
        return sqrt(max(ar.m^2, br.p^2) + max(az.m^2, bz.p^2))
    end
    
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

Internally, calls [𝒢_weno](@ref) on all computational cells.
"""
function 𝒢_weno_all(ϕ, dom::Domain)
    return reshape([𝒢_weno(ϕ, i, j, dom) for i in 1:dom.nr, j in 1:dom.nz], dom.nr, dom.nz)
end