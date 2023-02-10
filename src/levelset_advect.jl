export advect_ϕ, advect_ϕ!
"""
Upwind difference scheme, in the spirit of Godunov's scheme.
At the boundaries, if "upwind" goes outside domain, clamp gradients to 0.
"""
function calc_∇ϕ_2(ϕ, dom::Domain)
    dx1 = dom.dr1
    dy1 = dom.dz1
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    ∇ϕx = similar(ϕ)
    ∇ϕy = similar(ϕ)
    # outside_B = 1
    
    for j in 1:ny, i in 1:nx
        # At boundaries, if upwind is outside, clamp gradient to 0
        if i == 1  # Left edge
            pcell = ϕ[i,j]
            ∇ϕx[i,j] = flipsign(min(flipsign((ϕ[i+1,j] - ϕ[i,j] )*dx1, pcell), 0),pcell) 
        elseif i == nx # Right edge
            pcell = ϕ[i,j]
            ∇ϕx[i,j] = flipsign(max(flipsign((ϕ[i,j] - ϕ[i-1,j] )*dx1, pcell), 0),pcell) 
        else # Bulk
            pcell = flipsign(ϕ[i,j], ϕ[i,j])
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            if wcell < pcell && ecell < pcell # Conflux of two areas
                ∇ϕx[i,j] = (ϕ[i+1,j] - ϕ[i-1,j])*dx1*0.5 # Central difference
            elseif wcell < pcell
                ∇ϕx[i,j] = (ϕ[i,j] - ϕ[i-1,j] )*dx1   #Upwind to the left
            elseif ecell < pcell
                ∇ϕx[i,j] = (ϕ[i+1,j] - ϕ[i,j])*dx1 # Upwind to the right
            else # No upwind direction: clamp to 0
                ∇ϕx[i,j] = 0 # Upwind to the right
                # ∇ϕx[i,j] = (ϕ[i+1,j] - ϕ[i,j])*dx1 # Rightward difference
                
            end
        end
                
        if j == 1 # Bottom edge
            pcell = ϕ[i,j]
            ∇ϕy[i,j] = flipsign(min(flipsign((ϕ[i,j+1] - ϕ[i,j] )*dy1, pcell), 0),pcell) 
        elseif j == ny # Top edge
            pcell = ϕ[i,j]
            ∇ϕy[i,j] = flipsign(max(flipsign((ϕ[i,j] - ϕ[i,j-1] )*dy1, pcell), 0),pcell) 
        else # Bulk
            pcell = flipsign(ϕ[i,j], ϕ[i,j])
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell < pcell && ncell < pcell # Conflux of two areas
                ∇ϕy[i,j] = (ϕ[i,j+1] - ϕ[i,j-1])*dy1*0.5 # Centered difference
            elseif scell < pcell  # && scell < outside_B
                ∇ϕy[i,j] = (ϕ[i,j]- ϕ[i,j-1])*dy1 # Upwind downward
            elseif ncell < pcell # && ncell < outside_B
                ∇ϕy[i,j] = (ϕ[i,j+1] - ϕ[i,j])*dy1 # Upwind upward
            else # No upwind direction: clamp to 0
                ∇ϕy[i,j] = 0  
            end
        end
    end
    return ∇ϕx, ∇ϕy
end

"""
Upwind difference scheme for N ⋅ ∇ F , for each F ∈ v. Same finite differences for both N=∇ϕ and F.
At the boundaries, if no upwind direction available, clamp derivatives to 0.
"""
function calc_Nd∇v(ϕ, v, dom::Domain; debug=false)
    dx1 = dom.dr1
    dy1 = dom.dz1
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    dϕx = similar(ϕ)
    dϕy = similar(ϕ)
    dvx = similar(v)
    dvy = similar(v)
    for j in 1:ny, i in 1:nx
        # s = sign(ϕ[i,j]) # Use to ensure wind direction points away from contour
        pcell = flipsign(ϕ[i,j], ϕ[i,j])
        # At boundaries, if upwind not possible, clamp to 0
        if i == 1  # Left edge
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            if ecell > pcell# + 0.5dx
                dϕx[i,j]    = 0 # Clamp to 0 if  DOWNWIND
                dvx[i,j,:] .= 0 # Clamp to 0 if  DOWNWIND
            else
                dϕx[i,j]   =   (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                dvx[i,j,:] =@. (v[i+1,j,:] - v[i,j,:] )*dx1 
            end
        elseif i == nx # Right edge
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            if wcell > pcell# + 0.5dx
                dϕx[i,j]    = 0 # Clamp to 0 if  DOWNWIND
                dvx[i,j,:] .= 0 # Clamp to 0 if  DOWNWIND
            else
                dϕx[i,j]   =   (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                dvx[i,j,:] =@. (v[i,j,:] - v[i-1,j,:] )*dx1 
            end
        else # Bulk
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            wup = wcell < pcell
            eup = ecell < pcell
            if wup && eup
                dϕx[i,j]   =   (ϕ[i+1,j]   - ϕ[i-1,j]  )*0.5*dx1 # Second order central
                dvx[i,j,:] =@. (v[i+1,j,:] - v[i-1,j,:])*0.5*dx1 # Second order central
            elseif wup
                dϕx[i,j]   =   (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                dvx[i,j,:] =@. (v[i,j,:] - v[i-1,j,:] )*dx1 
            elseif eup 
                dϕx[i,j]   =   (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                dvx[i,j,:] =@. (v[i+1,j,:] - v[i,j,:] )*dx1 
            else # No upwind direction: clamp to 0.
                dϕx[i,j]    = 0
                dvx[i,j,:] .= 0
            end
        end
                
        if j == 1 # Bottom edge
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            if ncell > pcell# + 0.5dy
                dϕy[i,j]    = 0 # Clamp to 0 if would be DOWNWIND
                dvy[i,j,:] .= 0 # Clamp to 0 if would be DOWNWIND
            else
                dϕy[i,j]   =   (ϕ[i,j+1]   - ϕ[i,j]  )*dy1
                dvy[i,j,:] =@. (v[i,j+1,:] - v[i,j,:])*dy1
            end
        elseif j==ny # Top edge
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell > pcell# + 0.5dy
                dϕy[i,j]    = 0 # Clamp to 0 if would be DOWNWIND
                dvy[i,j,:] .= 0 # Clamp to 0 if would be DOWNWIND
            else
                dϕy[i,j]   =   (ϕ[i,j]   - ϕ[i,j-1]  )*dy1
                dvy[i,j,:] =@. (v[i,j,:] - v[i,j-1,:])*dy1
            end
        else # Bulk
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell < pcell && ncell < pcell 
                dϕy[i,j]   = (ϕ[i,j+1]   - ϕ[i,j-1  ])*0.5*dy1 # Second order central
                dvy[i,j,:] = (v[i,j+1,:] - v[i,j-1,:])*0.5*dy1 # Second order central
            elseif scell < pcell 
                dϕy[i,j]   =   (ϕ[i,j]   - ϕ[i,j-1]  )*dy1 # Upwind downward
                dvy[i,j,:] =@. (v[i,j,:] - v[i,j-1,:])*dy1 # Upwind upward
            elseif ncell < pcell
                dϕy[i,j]   =   (ϕ[i,j+1]   - ϕ[i,j]  )*dy1 # Upwind upward
                dvy[i,j,:] =@. (v[i,j+1,:] - v[i,j,:])*dy1 # Upwind upward
            else # No upwind direction: clamp to 0
                dϕy[i,j]    = 0
                dvy[i,j,:] .= 0
            end
        end
    end
    Nd∇v = @. dϕx * dvx + dϕy * dvy
    return Nd∇v
end

"""
    calc_Nd∇v!(cache, ϕ, v, dom::Domain; debug=false)

    Compute ∇ϕ ⋅ ∇F and store results in `cache`.

Upwind difference scheme for N ⋅ ∇ F , for each F ∈ v. Same finite differences for both N=∇ϕ and F.
At the boundaries, if no upwind direction available, clamp derivatives to 0.
"""
function calc_Nd∇v!(cache, ϕ, v, dom::Domain; debug=false)
    dx1 = dom.dr1
    dy1 = dom.dz1
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    nvec = size(v, 3)
    if size(cache) != size(v)
        @error "ArgumentError: improper cache size for storing Nd∇v"
    end
    # dϕx = similar(ϕ)
    # dϕy = similar(ϕ)
    # dvx = similar(v)
    # dvy = similar(v)
    for j in 1:ny, i in 1:nx
        # s = sign(ϕ[i,j]) # Use to ensure wind direction points away from contour
        pcell = flipsign(ϕ[i,j], ϕ[i,j])
        # At boundaries, if upwind not possible, clamp to 0
        if i == 1  # Left edge
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            if ecell > pcell# + 0.5dx
                dϕx = 0.0 # Clamp to 0 if  DOWNWIND
                dvx = fill(0.0, nvec) # Clamp to 0 if  DOWNWIND
            else
                dϕx =   (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                dvx =@. (v[i+1,j,:] - v[i,j,:] )*dx1 
            end
        elseif i == nx # Right edge
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            if wcell > pcell# + 0.5dx
                dϕx = 0.0 # Clamp to 0 if  DOWNWIND
                dvx = fill(0.0, nvec) # Clamp to 0 if  DOWNWIND
            else
                dϕx =   (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                dvx =@. (v[i,j,:] - v[i-1,j,:] )*dx1 
            end
        else # Bulk
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            wup = wcell < pcell
            eup = ecell < pcell
            if wup && eup
                dϕx =   (ϕ[i+1,j]   - ϕ[i-1,j]  )*0.5*dx1 # Second order central
                dvx =@. (v[i+1,j,:] - v[i-1,j,:])*0.5*dx1 # Second order central
            elseif wup
                dϕx =   (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                dvx =@. (v[i,j,:] - v[i-1,j,:] )*dx1 
            elseif eup 
                dϕx =   (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                dvx =@. (v[i+1,j,:] - v[i,j,:] )*dx1 
            else # No upwind direction: clamp to 0.
                dϕx = 0.0
                dvx = fill(0.0, nvec)
            end
        end
                
        if j == 1 # Bottom edge
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            if ncell > pcell# + 0.5dy
                dϕy = 0.0 # Clamp to 0 if would be DOWNWIND
                dvy = fill(0.0, nvec) # Clamp to 0 if would be DOWNWIND
            else
                dϕy =   (ϕ[i,j+1]   - ϕ[i,j]  )*dy1
                dvy =@. (v[i,j+1,:] - v[i,j,:])*dy1
            end
        elseif j==ny # Top edge
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell > pcell# + 0.5dy
                dϕy = 0.0 # Clamp to 0 if would be DOWNWIND
                dvy = fill(0.0, nvec) # Clamp to 0 if would be DOWNWIND
            else
                dϕy =   (ϕ[i,j]   - ϕ[i,j-1]  )*dy1
                dvy =@. (v[i,j,:] - v[i,j-1,:])*dy1
            end
        else # Bulk
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell < pcell && ncell < pcell 
                dϕy = (ϕ[i,j+1]   - ϕ[i,j-1  ])*0.5*dy1 # Second order central
                dvy = (v[i,j+1,:] - v[i,j-1,:])*0.5*dy1 # Second order central
            elseif scell < pcell 
                dϕy =   (ϕ[i,j]   - ϕ[i,j-1]  )*dy1 # Upwind downward
                dvy =@. (v[i,j,:] - v[i,j-1,:])*dy1 # Upwind upward
            elseif ncell < pcell
                dϕy =   (ϕ[i,j+1]   - ϕ[i,j]  )*dy1 # Upwind upward
                dvy =@. (v[i,j+1,:] - v[i,j,:])*dy1 # Upwind upward
            else # No upwind direction: clamp to 0
                dϕy = 0.0
                dvy = fill(0.0, nvec)
            end
        end
        @. cache[i,j,:] =  dϕx * dvx + dϕy * dvy
    end
    # Nd∇v = @. dϕx * dvx + dϕy * dvy
    # return Nd∇v
    nothing
end

"""
Takes vecfunc, a vector function which takes (ir, iz) and returns desired vector
"""
function vector_extrap_from_front(phi, Bf, vec_func, dom::Domain, dt=1.0, guess=nothing)
    # nr, nz = size(phi)
    nr, nz = dom.nr, dom.nz
    B = findall(Bf)
    front_cells = findall(identify_Γ(phi, dom) .& (phi.>= 0))
    nvec = size(vec_func(Tuple(front_cells[1])...), 1)
    if guess === nothing
        v0 = fill(0.0, nr, nz, nvec)
    else
        v0 = copy(guess)
    end
    # Include front cell values in initial condition, then leave them alone
    # maxval = 0
    for cell in front_cells
        val = vec_func(Tuple(cell)...)
        v0[cell,:] .= val
        # maxval = max(maxval, abs(val))
    end
    # scaled = sqrt(1/maxval)

    # ΩnB = findall(fill(true, nr, nz) .⊻ Bf)
    # v0[ΩnB] .= 0.0
    
    cached = copy(v0)
    # ndvcache = copy(v0)
    function sub_rhs(du, u, p, t)
        cached[B,:] .= reshape(u, :, nvec)
        # cached[B,:] .= u
        # F = reshape(u, nr, nz)

        # Thought an in-place might be faster, but had *more* allocations.
        # calc_Nd∇v!(ndvcache, phi, cached, dom;debug=false) #*scaled
        # ∂tv = - sign.(phi) .* ndvcache
        ∂tv = - sign.(phi) .* calc_Nd∇v( phi, cached, dom;debug=false) #*scaled

        for cell in front_cells
            ∂tv[cell,:] .= 0 # Don't mess with cells on the front
        end
        du .= reshape(∂tv[B,:], :)
        # return du
    end
    
    u0 = reshape(v0[B,:], :)
    CFL = 0.8
    tspan = (0.0, dt/CFL)
    # vnr, vnz = calc_∇ϕ_2(phi, dr, dz) # Not sure this is necessary.
    # subdt = CFL * minimum(@. dr / abs(vnr) + dz / abs(vnz))
    prob = ODEProblem(sub_rhs, u0, tspan)
    # sol = solve(prob, BS3(), dt=subdt; callback=TerminateSteadyState(1e-4, 1e-4))
    sol = solve(prob, BS3(); callback=TerminateSteadyState(1e-4, 1e-4))
    ret = cached
    ret[B,:] .= reshape(sol[end],:,nvec)
    return ret
end

    
"""
Given total speed v0, field ϕ and location (ir, iz), compute normal vector (into Ω⁻) times velocity v0.
"""
function compute_frontvel_1(v0, ϕ, i, j, dom::Domain)
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    dr1 = dom.dr1
    dz1 = dom.dz1
    pcell = flipsign(ϕ[i,j], ϕ[i,j])
    if pcell > 2dom.dr || pcell > 2dom.dz 
        @warn "Computed front velocity away from front: i=$i, j=$j, ϕ=$(ϕ[i,j])"
    end
    # At boundaries, if upwind not possible, clamp to 0
        # At boundaries, if upwind is outside, clamp gradient to 0
    if i == 1  # Left edge
        pcell = ϕ[i,j]
        ∇ϕx = flipsign(min(flipsign((ϕ[i+1,j] - ϕ[i,j] )*dr1, pcell), 0),pcell) 
    elseif i == nx # Right edge
        pcell = ϕ[i,j]
        ∇ϕx = flipsign(max(flipsign((ϕ[i,j] - ϕ[i-1,j] )*dr1, pcell), 0),pcell) 
    else # Bulk
        pcell = flipsign(ϕ[i,j], ϕ[i,j])
        ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
        wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
        if wcell < pcell && ecell < pcell # Conflux of two areas
            ∇ϕx = (ϕ[i+1,j] - ϕ[i-1,j])*dr1*0.5 # Central difference
        elseif wcell < pcell
            ∇ϕx = (ϕ[i,j] - ϕ[i-1,j] )*dr1   #Upwind to the left
        elseif ecell < pcell
            ∇ϕx = (ϕ[i+1,j] - ϕ[i,j])*dr1 # Upwind to the right
        else # No upwind direction: clamp to 0
            ∇ϕx = 0 
            
        end
    end
            
    if j == 1 # Bottom edge
        pcell = ϕ[i,j]
        ∇ϕy = flipsign(min(flipsign((ϕ[i,j+1] - pcell )*dz1, pcell), 0),pcell) 
    elseif j == ny # Top edge
        pcell = ϕ[i,j]
        ∇ϕy = flipsign(max(flipsign((pcell - ϕ[i,j-1] )*dz1, pcell), 0),pcell) 
    else # Bulk
        pcell = flipsign(ϕ[i,j], ϕ[i,j])
        ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
        scell = flipsign(ϕ[i,j-1],ϕ[i,j])
        if scell < pcell && ncell < pcell # Conflux of two areas
            ∇ϕy = (ϕ[i,j+1] - ϕ[i,j-1])*dz1*0.5 # Centered difference
        elseif scell < pcell  # && scell < outside_B
            ∇ϕy = (ϕ[i,j]- ϕ[i,j-1])*dz1 # Upwind downward
        elseif ncell < pcell # && ncell < outside_B
            ∇ϕy = (ϕ[i,j+1] - ϕ[i,j])*dz1 # Upwind upward
        else # No upwind direction: clamp to 0
            ∇ϕy = 0  
        end
    end
    return -v0*∇ϕx,  -v0*∇ϕy
end

"""
Upwind difference scheme for V⋅∇ϕ .
Vx = Vf[:,:,1], Vy = Vf[:,:,2] assumed.
Vx and Vy used to determine wind direction.
At the boundaries, move downwind where necessary, so there is at least an approximation for the gradient.
There are major performance gains to be had here--both in allocation and parallelization. TODO
"""
function calc_Vd∇ϕ(ϕ, Vf, dom ; outside_B = 1)
    dx1 = dom.dr1
    dy1 = dom.dz1
    # i12 = 1/12
    # i23 = 2/3
    Vx = Vf[:,:,1]
    Vy = Vf[:,:,2]
    # nx, ny = size(ϕ)
    nx = dom.nr
    ny = dom.nz
    ∇ϕx = similar(ϕ)
    ∇ϕy = similar(ϕ)
    Vd∇ϕ = similar(ϕ)

    for i in 1:nx, j in 1:ny
        pcell = ϕ[i,j]
        px = Vx[i,j] 
        py = Vy[i,j]
        outside_B = 1
        # Use velocity to dictate upwind; don't allow differences outside the computational band
        if px > 0
            if i == 1 
                if abs(ϕ[i+1,j]) == outside_B
                    ∇ϕx[i,j] = 0 
                else
                    ∇ϕx[i,j] = (ϕ[i+1,j] - pcell)*dx1 # DOWNWIND
                end
            else
                ∇ϕx[i,j] = (pcell - ϕ[i-1,j] )*dx1
            end
        elseif px < 0
            if i == nx
                # ∇ϕx[i,j] = (pcell - ϕ[i-1,j] )*dx1 # DOWNWIND
                # ∇ϕx[i,j] = 0 # DOWNWIND
                # ∇ϕy[i,j] = (abs(ϕ[i-1,j]) == outside_B ? 0 : (pcell - ϕ[i-1,j] )*dx1) # DOWNWIND
                if abs(ϕ[i-1,j]) == outside_B
                    ∇ϕx[i,j] = 0 
                    # println("clamped")
                else
                    ∇ϕx[i,j] = (pcell - ϕ[i-1,j] )*dx1 # DOWNWIND
                end
            else
                ∇ϕx[i,j] = (ϕ[i+1,j] - pcell)*dx1
            end
        else # Vx = 0: Vx *∇ϕx = 0
            ∇ϕx[i,j] = 0
        end
        if py > 0
            if j == 1 
                # ∇ϕy[i,j] = (ϕ[i,j+1] - pcell)*dy1 # DOWNWIND
                # ∇ϕy[i,j] = 0 # DOWNWIND
                # ∇ϕy[i,j] = (abs(ϕ[i,j+1]) == outside_B ? 0 : (ϕ[i,j+1] - pcell)*dy1) # DOWNWIND
                if abs(ϕ[i,j+1]) == outside_B
                    ∇ϕy[i,j] =  0 
                    # println("clamped")
                else
                    ∇ϕy[i,j] = (ϕ[i,j+1] - pcell)*dy1 # DOWNWIND
                end
            else
                ∇ϕy[i,j] = (pcell - ϕ[i,j-1])*dy1
            end
        elseif py < 0
            if j == ny
                # ∇ϕy[i,j] = (pcell - ϕ[i,j-1])*dy1 # DOWNWIND
                # ∇ϕy[i,j] = 0 # DOWNWIND
                # ∇ϕy[i,j] = (abs(ϕ[i,j-1]) == outside_B ? 0 : (pcell - ϕ[i,j-1] )*dy1) # DOWNWIND
                if abs(ϕ[i,j-1]) == outside_B
                    ∇ϕy[i,j] =  0 
                    # println("clamped")
                else
                    ∇ϕy[i,j] = (pcell - ϕ[i,j-1])*dy1 # DOWNWIND
                end
            else
                ∇ϕy[i,j] = (ϕ[i,j+1] - pcell)*dy1
            end
        else # Vy = 0: no gradient contribution
            ∇ϕy[i,j] = 0
        end
        
    end
    @. Vd∇ϕ = Vx * ∇ϕx + Vy * ∇ϕy
    return Vd∇ϕ
end

function advect_ϕ(ϕ, Vf, dom::Domain, dt)
    # nr, nz = size(ϕ)
    @unpack nr, nz = dom
    # nr = dom.nr
    # nz = dom.nz
    function sub_rhs(u, p, t)
        ϕ = reshape(u, nr, nz)
        ∂tϕ = -calc_Vd∇ϕ(ϕ, Vf, dom)
        du = reshape(∂tϕ, :)
        return du
    end
    u0 = reshape(ϕ, :)
    tspan = (0.0, dt)
    CFL = 0.8
    Vr = Vf[:,:,1]
    Vz = Vf[:,:,1]
    if maximum(abs.(Vr)) == 0
        subdt = CFL * minimum(@. dom.dz / abs(Vz))
    elseif maximum(abs.(Vz)) == 0
        subdt = CFL * minimum(@. dom.dr / abs(Vr))
    else
        subdt = CFL * minimum(@. dom.dr / abs(Vr) + dom.dz / abs(Vz))
    end
    prob = ODEProblem(sub_rhs, u0, tspan)
    sol = solve(prob, SSPRK43(), dt=subdt)
    # display(sol)
    # display(sol[end] .- u0)
    return reshape(sol[end], nr, nz)
end
function advect_ϕ!(ϕ, Vf, dom::Domain, dt)
    ϕ .= advect_ϕ(ϕ, Vf, dom, dt)
    return nothing
end