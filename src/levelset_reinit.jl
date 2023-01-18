function identify_Γ(ϕ)
    locs = similar(ϕ, Bool)
    locs .= false
    sg = sign.(ϕ)
    xshift = sg[1:end-1,:] .* sg[2:end,:]
    yshift = sg[:,1:end-1] .* sg[:,2:end]
    nx, ny = size(ϕ)
    for i in 1:nx-1, j in 1:ny
        # if locs[i,j]
        #     continue
        # end
        if xshift[i,j] <= 0
            locs[i,j] = locs[i+1,j] = true
        end
    end
    for i in 1:nx, j in 1:ny-1
        # if locs[i,j]
        #     continue
        # end
        if yshift[i,j] <= 0
            locs[i,j] = locs[i, j+1] = true
        end
    end
    return locs
end

Γ_cells(ϕ) = findall(identify_Γ(ϕ))

"Not used explicitly here, just useful for debugging."
function calc_curvature(ϕ, dx, dy)
    dx2 = 1/dx^2
    dy2 = 1/dy^2
    nx, ny = size(ϕ)
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
Takes full level set field ϕ, list of front cells Γ, dx, and dy.
Computes curvature (or at least something proportional to it) at all locations Γ, then compares against sign of ϕ to assign to R or C
"""
function identify_regions_RC(ϕ, Γ, dx, dy)
    dx2 = 1/dx^2
    dy2 = 1/dy^2
    nx, ny = size(ϕ)
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

# function plot_RC(RC, nx, ny)
#     R, C = RC
#     arr = fill(0, nx, ny)
#     for c in R
#         arr[c] += 1
#     end
#     for c in C
#         arr[c] -= 1
#     end
#     heat(arr)
# end
function plot_RC(ϕ)
    R, C = identify_regions_RC(ϕ, Γ_cells(ϕ), dr, dz)
    Rr = [rgrid[Tuple(c)[1]] for c in R]
    Rz = [zgrid[Tuple(c)[2]] for c in R]
    Cr = [rgrid[Tuple(c)[1]] for c in C]
    Cz = [zgrid[Tuple(c)[2]] for c in C]
    scatter!(Rr, Rz, c=:black)
    scatter!(Cr, Cz, c=:white)
    # plot_RC(RC, nr, nz)
end

"""
Takes a field of bools identifying Γ, bandwidth in x cells, and bandwidth in y cells
Returns a field of bools identifying B
"""
function identify_B(Γ_field::Matrix{Bool}, bwx, bwy)
    nx, ny = size(Γ_field)
    B = fill(false, nx, ny)
    Γ = findall(Γ_field)
    for c in Γ
        ix, iy = Tuple(c)
        xgrab = range(max(1,ix-bwx), min(nx, ix+bwx))
        ygrab = range(max(1,iy-bwy), min(ny, iy+bwy))
        B[xgrab, iy] .= true
        B[ix, ygrab] .= true
    end
    return B
end
"""
Takes a field of bools identifying Γ, bandwidth in x cells, and bandwidth in y cells
Returns a field of bools identifying B
"""
function identify_B(ϕ::Matrix{Float64}, bwx, bwy)
    Γ_field = identify_Γ(ϕ)
    nx, ny = size(Γ_field)
    B = fill(false, nx, ny)
    Γ = findall(Γ_field)
    for c in Γ
        ix, iy = Tuple(c)
        xgrab = range(max(1,ix-bwx), min(nx, ix+bwx))
        ygrab = range(max(1,iy-bwy), min(ny, iy+bwy))
        B[xgrab, iy] .= true
        B[ix, ygrab] .= true
    end
    return B
end

"""
Uses global rgrid, dr
"""
function calc_dϕdr_sdf(ϕ, Γf, i, j, nr, nz)
    # 
    if i == nr
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
    den = max(rgrid[ip] - rgrid[im], .001*dr)
    return num/den
end

"""
Uses global zgrid, dz
"""
function calc_dϕdz_sdf(ϕ, Γf, i, j, nr, nz)
    # Fancy conditions for near coalescence: ignored for now
    # TODO: fill these out for real. If A && B is true, jp = j, likewise for jm = j or something
    A = true
    B = true
    if j == nz
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
    den = max(zgrid[jp] - zgrid[jm], .001*dz)
    return num/den
end

function calc_dij_R!(d, ϕ, Γf, R)
    nr, nz = size(ϕ)
    for c in R
        if ϕ[c]==0
            d[c] = 0
            continue
        end
        i, j = Tuple(c)
        den = hypot(calc_dϕdr_sdf(ϕ, Γf, i, j, nr, nz), calc_dϕdz_sdf(ϕ, Γf, i, j, nr, nz))
        d[c] = ϕ[c] / den
    end
end
function calc_dij!(d, ϕ, Γf, R)
    nr, nz = size(ϕ)
    for c in R
        if ϕ[c]==0
            d[c] = 0
            continue
        end
        i, j = Tuple(c)
        den = hypot(calc_dϕdr_sdf(ϕ, Γf, i, j, nr, nz), calc_dϕdz_sdf(ϕ, Γf, i, j, nr, nz))
        d[c] = ϕ[c] / den
    end
end
function calc_dij_C!(d, ϕ, Γf, C)
    nr, nz = size(ϕ)
    for c in C
        if(ϕ[c]==0)
            d[c] = 0
            continue
        end
        i, j = Tuple(c)
        neighbors = Vector{Tuple}()
        if i == 1
            push!(neighbors, (i+1,j))
        elseif i == nr
            push!(neighbors, (i-1,j))
        else
            push!(neighbors, (i+1,j))
            push!(neighbors, (i-1,j))
        end
        if j == 1
            push!(neighbors, (i,j+1))
        elseif j == nz
            push!(neighbors, (i,j-1))
        else
            push!(neighbors, (i,j+1))
            push!(neighbors, (i,j-1))
        end

        Sij = [nb for nb in neighbors if ϕ[nb...]*ϕ[c] < 0]
        if length(Sij) > 0
            num = sum([d[nb...] for nb in Sij])
            den = sum([ϕ[nb...] for nb in Sij])
        else
            # Happens because a cell is exactly 0, so Γ is three cells wide.
            # Identify the 0 neighbor, set distance to neighbor
            # println("Watch for this case! ")
            println("Length of Sij is $(length(Sij))")
            Sij = [nb for nb in neighbors if ϕ[nb...] == 0]
            if (i + 1,j) ∈ Sij || (i-1,j) ∈ Sij
                mindx = dr
            end
            if (i,j+1) ∈ Sij || (i,j-1) ∈ Sij
                mindx = dz
            end
            d[c] = sign(ϕ[c]) * mindx
            continue
            # num = den = 1 #?
            
            # println(ϕ[c])
        end
        d[c] =  ϕ[c] * num / den
    end
end
function calc_dtldij(d, ϕ, Γf, cell)
    nr, nz = size(ϕ)
    if(ϕ[cell]==0)
        return 0
    end
    i, j = Tuple(cell)
    neighbors = Vector{Tuple}()
    if i == 1
        push!(neighbors, (i+1,j))
    elseif i == nr
        push!(neighbors, (i-1,j))
    else
        push!(neighbors, (i+1,j))
        push!(neighbors, (i-1,j))
    end
    if j == 1
        push!(neighbors, (i,j+1))
    elseif j == nz
        push!(neighbors, (i,j-1))
    else
        push!(neighbors, (i,j+1))
        push!(neighbors, (i,j-1))
    end

    Sij = [nb for nb in neighbors if ϕ[nb...]*ϕ[cell] < 0]
    if length(Sij) > 0
        num = sum([d[nb...] for nb in Sij])
        den = sum([ϕ[nb...] for nb in Sij])
    else
        # Happens because a cell is exactly 0, so Γ is three cells wide.
        # Identify the 0 neighbor, set distance to neighbor
        # println("Watch for this case! ")
        println("Length of Sij is $(length(Sij))")
        Sij = [nb for nb in neighbors if ϕ[nb...] == 0]
        if (i + 1,j) ∈ Sij || (i-1,j) ∈ Sij
            mindx = dr
        end
        if (i,j+1) ∈ Sij || (i,j-1) ∈ Sij
            mindx = dz
        end
        return sign(ϕ[cell]) * mindx
        # num = den = 1 #?
        
        # println(ϕ[c])
    end
    return ϕ[c] * num / den
end
"Uses global dr, dz"
function update_ϕ_in_Γ!(ϕl)
    Γfl = identify_Γ(ϕl)
    Γl = findall(Γfl)
    RCl = identify_regions_RC(ϕl, Γl, dr, dz)
    nr, nz = size(ϕl)
    dl = fill(0.0, nr, nz)
    calc_dij_R!(dl, ϕl, Γfl, Γl)
    # dl2 = copy(dl)
    calc_dij_C!(dl, ϕl, Γfl, RCl[2])
    calc_dij_C!(dl, ϕl, Γfl, RCl[2])
    # calc_dij!(dl, ϕl, Γfl, Γl)
    # dtld = calc_dtldij(dl, ϕl, Γfl, RCl[2])
    for c in Γl
        ϕl[c] = dl[c]
        # if c ∈ RCl[2]
        #     @show(calc_dtldij(dl2,ϕl,Γfl,c) - dl[c])
        # end
        #     ϕl[c] = calc_dtldij(dl, ϕl, Γfl, c)
        # else
        #     ϕl[c] = dl[c]
        # end
    end
    # for i in 1:nr, j in 1:nz
    #     if Γf[i,j]
    #         ϕ[i,j] = d[i,j]
    #     end
    # end
    # arr = fill(0.0, nr, nz)
    # arr[Γ] .= du
    # display(heat(d))
    # return dl
end

struct LD{T} # LD short for Little Difference
    p::T
    m::T
end
LD(x) = LD(max(x, 0), min(x, 0))
"""
Godunov's scheme for discretizing the norm of the gradient of ϕ.
Uses global nr, nz, dr, dz
"""
function 𝒢(ϕ, i, j) # p. 6830 of Hartmann, 10th page of PDF
    dr1 = 1/dr
    dz1 = 1/dz
    pcell = ϕ[i,j]
    if i == 1
        a = LD(0)
        b = LD((ϕ[i+1,j] - ϕ[i,j]) * dr1)
    elseif i == nr
        a = LD((ϕ[i,j] - ϕ[i-1,j]) * dr1)
        b = LD(0)
    else
        a = LD((ϕ[i,j] - ϕ[i-1,j]) * dr1)
        b = LD((ϕ[i+1,j] - ϕ[i,j]) * dr1)
    end
    if j == 1
        c = LD(0)
        d = LD((ϕ[i,j+1] - ϕ[i,j]) * dz1)
    elseif j == nz
        c = LD((ϕ[i,j] - ϕ[i,j-1]) * dz1)
        d = LD(0)
    else
        c = LD((ϕ[i,j] - ϕ[i,j-1]) * dz1)
        d = LD((ϕ[i,j+1] - ϕ[i,j]) * dz1)
    end
    if ϕ[i,j] >= 0
    # if ϕ[i,j] < 0
        return sqrt(max(a.p^2, b.m^2) + max(c.p^2, d.m^2))
    else
        return sqrt(max(a.m^2, b.p^2) + max(c.m^2, d.p^2))
    end
end

function 𝒢_all(ϕ)
    return reshape([𝒢(ϕ, i, j) for i in 1:nr, j in 1:nz], nr, nz)
end


function reinitialize_ϕ!(ϕ_mat, tf=1.0; alg=BS3(), bwr=5, bwz=4)

    Γf = identify_Γ(ϕ_mat)
    Γ = findall(Γf)
    # bwr = 5
    # bwz = 4
    Bf = identify_B(Γf, bwr, bwz)
    BnΓ = findall(Bf .⊻ Γf)
    ΩnB = findall(fill(true, nr, nz) .⊻ Bf)

    update_ϕ_in_Γ!(ϕ_mat)

    nx, ny = size(ϕ_mat)
    sarr = sign.(ϕ_mat)
    Γ = Γ_cells(ϕ_mat)

    
    # ϕ_ode = reshape(ϕ_mat, :)
    ϕ_ode = ϕ_mat[BnΓ]
    cached = copy(ϕ_mat)
    function sub_rhs(du, u, p, t) 
        cached[BnΓ] .= u
        # dϕ = sarr .* (1 .- 𝒢_all(cached))
        # return dϕ[BnΓ]
        # du = zeros(length(BnΓ))
        for (i, c) in enumerate(BnΓ)
            du[i] = sarr[c] * (1-𝒢(cached, Tuple(c)...))
        end
        return du


        # dϕ = sarr .* (1 .- 𝒢_all(reshape(u, nx, ny)))
        # dϕ[Γ] .= 0.0
        # dϕ[ΩnB] .= 0.0
        # return reshape(dϕ, :)
    end
    tspan = (0.0, tf)
    prob = ODEProblem(sub_rhs, ϕ_ode, tspan)
    # sol = solve(prob, Tsit5(), dt = 1.0; callback=TerminateSteadyState(1e-4, 1e-4))
    # sol = solve(prob, Euler(), dt = 1.0; callback=TerminateSteadyState(1e-4, 1e-4))
    sol = solve(prob, alg, dt = 1.0; callback=TerminateSteadyState(1e-4, 1e-4))
    # ϕ_sol = ϕ_mat
    # ϕ_sol[BnΓ] .= sol[end]
    ϕ_mat[BnΓ] .= sol[end]

    # ϕ_sol[ΩnB] .= sarr[ΩnB]
    ϕ_mat[ΩnB] .= sarr[ΩnB]

    # ϕ_sol
    nothing
    # ϕ_rep = reshape(sol[end], nx, ny)
end

function reinitialize_ϕ(ϕ, tf=1.0)
    ϕ1 = copy(ϕ)
    reinitialize_ϕ!(ϕ1, tf)
    ϕ1
end