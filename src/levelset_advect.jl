export advect_ϕ, advect_ϕ!
export vector_extrap_from_front
export fastmarch_v!, extrap_v_fastmarch, extrap_v_pde

# FUnctions exported for Documenter.jl sake
export calc_Vd∇ϕ
export calc_Nd∇v, calc_Nd∇v!
export calc_∇ϕ_1st, compute_frontvel_1

"""
    calc_∇ϕ_1st(ϕ, dom::Domain)

Compute ∇ϕ with a 1st-order "upwind" difference scheme ("upwind" meaning towards interface).

At present, this is not being used
At the boundaries, if "upwind" goes outside domain, clamp derivatives to 0.

"""
function calc_∇ϕ_1st(ϕ, dom::Domain)
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
    calc_Nd∇v(ϕ, v, dom::Domain; debug=false)

Compute ∇ϕ ⋅ ∇F for each F ∈ v with a first-order upwind scheme ("upwind" toward interface.).

Full allocating.
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
    # Nd∇v = @. dϕx * dvx + dϕy * dvy
    return @. dϕx * dvx + dϕy * dvy
end

"""
    calc_Nd∇v!(cache, ϕ, v, dom::Domain; debug=false)

Compute ∇ϕ ⋅ ∇F and store results in `cache`.

Stores all results in `cache`; tries to avoid allocation.
1st order Upwind difference scheme for N ⋅ ∇ F , for each F ∈ v. Same finite differences for both N=∇ϕ and F.
At the boundaries, if no upwind direction available, clamp derivatives to 0.
"""
function calc_Nd∇v!(cache, ϕ, v, dom::Domain)
    dx1 = dom.dr1
    dy1 = dom.dz1
    nvec = size(v, 3)
    if size(cache) != size(v)
        @error "ArgumentError: improper cache size for storing Nd∇v"
    end
    # dϕx = similar(ϕ)
    # dϕy = similar(ϕ)
    # dvx = similar(v)
    # dvy = similar(v)
    # Cache some variables for saving on memory
    dϕx = 0.0
    dϕy = 0.0
    dvx = fill(0.0, nvec)
    dvy = fill(0.0, nvec)
    for j in 1:dom.nz, i in 1:dom.nr

        # s = sign(ϕ[i,j]) # Use to ensure wind direction points away from contour
        pcell = flipsign(ϕ[i,j], ϕ[i,j])
        # At boundaries, if upwind not possible, clamp to 0
        if i == 1  # Left edge
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            if ecell > pcell# + 0.5dx
                   dϕx = 0.0 # Clamp to 0 if  DOWNWIND
                @. dvx = fill(0.0, nvec) # Clamp to 0 if  DOWNWIND
            else
                   dϕx = (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                @. dvx = (v[i+1,j,:] - v[i,j,:] )*dx1 
            end
        elseif i == dom.nr # Right edge
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            if wcell > pcell# + 0.5dx
                   dϕx = 0.0 # Clamp to 0 if  DOWNWIND
                @. dvx = fill(0.0, nvec) # Clamp to 0 if  DOWNWIND
            else
                   dϕx = (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                @. dvx = (v[i,j,:] - v[i-1,j,:] )*dx1 
            end
        else # Bulk
            ecell = flipsign(ϕ[i+1,j],ϕ[i,j])
            wcell = flipsign(ϕ[i-1,j],ϕ[i,j])
            wup = wcell < pcell
            eup = ecell < pcell
            if wup && eup
                   dϕx = (ϕ[i+1,j]   - ϕ[i-1,j]  )*0.5*dx1 # Second order central
                @. dvx = (v[i+1,j,:] - v[i-1,j,:])*0.5*dx1 # Second order central
            elseif wup
                   dϕx = (ϕ[i,j]   - ϕ[i-1,j]   )*dx1 
                @. dvx = (v[i,j,:] - v[i-1,j,:] )*dx1 
            elseif eup 
                   dϕx = (ϕ[i+1,j]   - ϕ[i,j]   )*dx1 
                @. dvx = (v[i+1,j,:] - v[i,j,:] )*dx1 
            else # No upwind direction: clamp to 0.
                   dϕx = 0.0
                @. dvx = fill(0.0, nvec)
            end
        end
                
        if j == 1 # Bottom edge
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            if ncell > pcell# + 0.5dy
                   dϕy = 0.0 # Clamp to 0 if would be DOWNWIND
                @. dvy = fill(0.0, nvec) # Clamp to 0 if would be DOWNWIND
            else
                   dϕy = (ϕ[i,j+1]   - ϕ[i,j]  )*dy1
                @. dvy = (v[i,j+1,:] - v[i,j,:])*dy1
            end
        elseif j==dom.nz # Top edge
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell > pcell# + 0.5dy
                   dϕy = 0.0 # Clamp to 0 if would be DOWNWIND
                @. dvy = fill(0.0, nvec) # Clamp to 0 if would be DOWNWIND
            else
                   dϕy = (ϕ[i,j]   - ϕ[i,j-1]  )*dy1
                @. dvy = (v[i,j,:] - v[i,j-1,:])*dy1
            end
        else # Bulk
            ncell = flipsign(ϕ[i,j+1],ϕ[i,j])
            scell = flipsign(ϕ[i,j-1],ϕ[i,j])
            if scell < pcell && ncell < pcell 
                   dϕy = (ϕ[i,j+1]   - ϕ[i,j-1  ])*0.5*dy1 # Second order central
                @. dvy = (v[i,j+1,:] - v[i,j-1,:])*0.5*dy1 # Second order central
            elseif scell < pcell 
                   dϕy = (ϕ[i,j]   - ϕ[i,j-1]  )*dy1 # Upwind downward
                @. dvy = (v[i,j,:] - v[i,j-1,:])*dy1 # Upwind upward
            elseif ncell < pcell
                   dϕy = (ϕ[i,j+1]   - ϕ[i,j]  )*dy1 # Upwind upward
                @. dvy = (v[i,j+1,:] - v[i,j,:])*dy1 # Upwind upward
            else # No upwind direction: clamp to 0
                   dϕy = 0.0
                @. dvy = fill(0.0, nvec)
            end
        end
        @. cache[i,j,:] =  dϕx * dvx + dϕy * dvy
    end
    # Nd∇v = @. dϕx * dvx + dϕy * dvy
    # return Nd∇v
    nothing
end

"""
    calc_Nd∇v_noalloc(ir, iz, ϕ, v, dom::Domain, cache)

Compute ∇ϕ ⋅ ∇F at `ir` and `iz`, using `cache` to store all intermediate values.

1st order Upwind difference scheme for N ⋅ ∇ F , for each F ∈ v. Same finite differences for both N=∇ϕ and F.
At the boundaries, if no upwind direction available, clamp derivatives to 0.
"""
function calc_Nd∇v_noalloc!(ndvcache, ϕ, v, dom::Domain, intermcache)
    dx1 = dom.dr1
    dy1 = dom.dz1
    nvec = size(v, 3)
    if size(intermcache, 1) < 7+2nvec
        @error "ArgumentError: too-small cache size for holding Nd∇v intermediates" 
    end
    # Cached variables for saving on memory
    pc, ec, wc, sc, nc = @view intermcache[1:5]
    dϕx = @view intermcache[6]
    dϕy = @view intermcache[7]
    dvx = @view intermcache[8:7+nvec]
    dvy = @view intermcache[8+nvec:7+2nvec]


    for iz in 1:dom.nz, ir in 1:dom.nr
        # s = sign(ϕ[ir,iz]) # Use to ensure wind direction points away from contour
        pc = flipsign(ϕ[ir,iz], ϕ[ir,iz])
        # At boundaries, if upwind not possible, clamp to 0
        if ir == 1  # Left edge
            ec = flipsign(ϕ[ir+1,iz],ϕ[ir,iz])
            if ec > pc# + 0.5dx
                dϕx .= 0.0 # Clamp to 0 if  DOWNWIND
                dvx .= 0.0 # Clamp to 0 if  DOWNWIND
            else
                dϕx .= (ϕ[ir+1,iz]   - ϕ[ir,iz]   )*dx1 
                dvx .= (v[ir+1,iz,:] - v[ir,iz,:] )*dx1 
            end
        elseif ir == dom.nr # Right edge
            wc = flipsign(ϕ[ir-1,iz],ϕ[ir,iz])
            if wc > pc# + 0.5dx
                dϕx .= 0.0 # Clamp to 0 if  DOWNWIND
                dvx .= 0.0 # Clamp to 0 if  DOWNWIND
            else
                dϕx .= (ϕ[ir,iz]   - ϕ[ir-1,iz]   )*dx1 
                dvx .= (v[ir,iz,:] - v[ir-1,iz,:] )*dx1 
            end
        else # Bulk
            ec = flipsign(ϕ[ir+1,iz],ϕ[ir,iz])
            wc = flipsign(ϕ[ir-1,iz],ϕ[ir,iz])
            if ec<pc && wc<pc
                dϕx .= (ϕ[ir+1,iz]   - ϕ[ir-1,iz]  )*0.5*dx1 # Second order central
                dvx .= (v[ir+1,iz,:] - v[ir-1,iz,:])*0.5*dx1 # Second order central
            elseif wc<pc
                dϕx .= (ϕ[ir,iz]   - ϕ[ir-1,iz]   )*dx1 
                dvx .= (v[ir,iz,:] - v[ir-1,iz,:] )*dx1 
            elseif ec<pc 
                dϕx .= (ϕ[ir+1,iz]   - ϕ[ir,iz]   )*dx1 
                dvx .= (v[ir+1,iz,:] - v[ir,iz,:] )*dx1 
            else # No upwind direction: clamp to 0.
                dϕx .= 0.0
                dvx .= 0.0
            end
        end
                
        if iz == 1 # Bottom edge
            nc = flipsign(ϕ[ir,iz+1],ϕ[ir,iz])
            if nc > pc# + 0.5dy
                dϕy .= 0.0 # Clamp to 0 if would be DOWNWIND
                dvy .= 0.0 # Clamp to 0 if would be DOWNWIND
            else
                dϕy .= (ϕ[ir,iz+1]   - ϕ[ir,iz]  )*dy1
                dvy .= (v[ir,iz+1,:] - v[ir,iz,:])*dy1
            end
        elseif iz==dom.nz # Top edge
            sc = flipsign(ϕ[ir,iz-1],ϕ[ir,iz])
            if sc > pc# + 0.5dy
                dϕy .= 0.0 # Clamp to 0 if would be DOWNWIND
                dvy .= 0.0 # Clamp to 0 if would be DOWNWIND
            else
                dϕy .= (ϕ[ir,iz]   - ϕ[ir,iz-1]  )*dy1
                dvy .= (v[ir,iz,:] - v[ir,iz-1,:])*dy1
            end
        else # Bulk
            nc = flipsign(ϕ[ir,iz+1],ϕ[ir,iz])
            sc = flipsign(ϕ[ir,iz-1],ϕ[ir,iz])
            if sc < pc && nc < pc 
                dϕy .= (ϕ[ir,iz+1]   - ϕ[ir,iz-1  ])*0.5*dy1 # Second order central
                dvy .= (v[ir,iz+1,:] - v[ir,iz-1,:])*0.5*dy1 # Second order central
            elseif sc < pc 
                dϕy .= (ϕ[ir,iz]   - ϕ[ir,iz-1]  )*dy1 # Upwind downward
                dvy .= (v[ir,iz,:] - v[ir,iz-1,:])*dy1 # Upwind upward
            elseif nc < pc
                dϕy .= (ϕ[ir,iz+1]   - ϕ[ir,iz]  )*dy1 # Upwind upward
                dvy .= (v[ir,iz+1,:] - v[ir,iz,:])*dy1 # Upwind upward
            else # No upwind direction: clamp to 0
                dϕy .= 0.0
                dvy .= 0.0
            end
        end
        # return dϕx * dvx + dϕy * dvy
        @. ndvcache[ir,iz,:] = dϕx * dvx + dϕy * dvy
    end

    # Nd∇v = @. dϕx * dvx + dϕy * dvy
    # return Nd∇v
    nothing
end

"""
Takes vecfunc, a vector function which takes (ir, iz) and returns desired vector
"""
function vector_extrap_from_front(phi, Bf, vec_func, dom::Domain, dt=1.0, guess=nothing)
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
        v0[cell,:] .= vec_func(Tuple(cell)...)
        # maxval = max(maxval, abs(val))
    end
    # scaled = sqrt(1/maxval)

    # ΩnB = findall(fill(true, nr, nz) .⊻ Bf)
    # v0[ΩnB] .= 0.0
    
    vcache = copy(v0)
    ndvcache = copy(v0)
    intermcache = fill(0.0, 7 + 2*2) # nvec = 2
    function sub_rhs(du, u, p, t)
        vcache[B,:] .= reshape(u, :, nvec)
        # vcache[B,:] .= u
        # F = reshape(u, nr, nz)

        # Thought an in-place might be faster, but had *more* allocations.
        calc_Nd∇v_noalloc!(ndvcache, phi, vcache, dom, intermcache) #*scaled
        ∂tv = - sign.(phi) .* ndvcache
        # ∂tv = - sign.(phi) .* calc_Nd∇v( phi, cached, dom;debug=false) #*scaled

        for cell in front_cells
            ∂tv[cell,:] .= 0 # Don't mess with cells on the front
        end
        du .= reshape(∂tv[B,:], :)
        # return du
    end
    
    u0 = reshape(v0[B,:], :)
    CFL = 0.8
    tspan = (0.0, dt/CFL)
    # vnr, vnz = calc_∇ϕ_1st(phi, dr, dz) # Not sure this is necessary.
    # subdt = CFL * minimum(@. dr / abs(vnr) + dz / abs(vnz))
    prob = ODEProblem(sub_rhs, u0, tspan)
    # sol = solve(prob, BS3(), dt=subdt; callback=TerminateSteadyState(1e-4, 1e-4))
    sol = solve(prob, BS3(); callback=TerminateSteadyState(1e-4, 1e-4))
    ret = vcache
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
    Vx = Vf[:,:,1]
    Vy = Vf[:,:,2]
    nx = dom.nr
    ny = dom.nz
    ∇ϕx = similar(ϕ)
    ∇ϕy = similar(ϕ)
    Vd∇ϕ = similar(ϕ)

    for i in 1:nx, j in 1:ny
        pcell = ϕ[i,j]
        px = Vx[i,j] 
        py = Vy[i,j]
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

"""
    advect_ϕ(ϕ, Vf, dom::Domain, dt; alg=SSPRK43())

Advect the level set field `ϕ` by velocity vield `Vf`, with time step `dt`.

`Vf` should have size `(nr, nz, 2)`. Time steps for the algorithm are taken using `alg`.
"""
function advect_ϕ(ϕ, Vf, dom::Domain, dt; alg=SSPRK43())
    @unpack nr, nz = dom
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
    sol = solve(prob, alg, dt=subdt)
    # display(sol)
    # display(sol[end] .- u0)
    return reshape(sol[end], nr, nz)
end
"""
    function advect_ϕ!(ϕ, Vf, dom::Domain, dt)

Thin wrapper on advect_ϕ--this is not necessarily better for performance.
"""
function advect_ϕ!(ϕ, Vf, dom::Domain, dt)
    ϕ = advect_ϕ(ϕ, Vf, dom, dt)
    return nothing
end


function ddx_fastmarch_2nd!(num, den, ϕ, vf, c, sdir, dersign, dx1)
    dϕdx = dersign*(-1.5ϕ[c] + 2ϕ[c+sdir] - 0.5ϕ[c+2sdir])*dx1
    px = -1.5*dersign*dx1 
    dv_x = dersign*(        2vf[c+sdir, :] .-0.5vf[c+2sdir,:]).*dx1

    @. num += dϕdx*dv_x
    den .+= dϕdx * px
end
function ddx_fastmarch_1st!(num, den, ϕ, vf, c, sdir, dersign, dx1)
    dϕdx = dersign*(-ϕ[c] + ϕ[c+sdir])*dx1
    px = -dersign*dx1 
    dv_x = dersign*(        vf[c+sdir, :]).*dx1

    @. num += dϕdx*dv_x
    den .+= dϕdx * px
end

"""
    fastmarch_v!(vf, acc, locs, ϕ, dom)

Mutate `vf` and `acc` to extrapolate `vf` from interface by fast marching.

Cells are calculated in order of increasing |ϕ|.
Uses matrix of bools `acc` (denoting accepted cells) to determine cells to use in extrapolation,
so if you determine one side first, will get used in the derivatives for the other side.
"""
function fastmarch_v!(vf, acc, locs::Vector{CartesianIndex{2}}, ϕ, dom::Domain)

    # Set up stencil checking tools

    usecell(acc, c) = checkbounds(Bool, acc, c) && acc[c]

    eci  = CartesianIndex( 1, 0)
    wci  = CartesianIndex(-1, 0)
    nci  = CartesianIndex( 0, 1)
    sci  = CartesianIndex( 0,-1)

    sort!(locs , by=(x->abs(ϕ[x])))

    for c in locs # This needs to be done in order.

        num = fill(0.0, 2)
        den = [0.0]

        # R direction
        if usecell(acc, c+eci) # Use east
            if usecell(acc, c+2eci) # Second order east
                # dϕdr = (-1.5ϕ[c] + 2ϕ[c+eci] - 0.5ϕ[c+e⁺ci])*dom.dr1
                # pr = -1.5dom.dr1 
                # dv_r = (        2vf[c+eci, :] .-0.5vf[c+e⁺ci,:]).*dom.dr1

                # @. num += dϕdr*dv_r
                # den += dϕdr * pr
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, eci, 1, dom.dr1)
                # println("East, 2nd order, c=$c")
            else
                # dϕdr = (-ϕ[c] + ϕ[c+eci])*dom.dr1
                # pr = -dom.dr1
                # dv_r =        vf[c+eci, :] .*dom.dr1

                # @. num += dϕdr*dv_r
                # den += dϕdr * pr
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, eci, 1, dom.dr1)
                # println("East, 1st order, c=$c")
            end
        elseif usecell(acc, c+wci)
            if usecell(acc, c+2wci) # Second order west
                # dϕdr = ( 1.5ϕ[c] - 2ϕ[c+wci] + 0.5ϕ[c+w⁻ci])*dom.dr1
                # pr =  1.5dom.dr1 
                # dv_r =       (-2vf[c+wci, :] .+0.5vf[c+w⁻ci,:]).*dom.dr1

                # @. num += dϕdr*dv_r
                # den += dϕdr * pr
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, wci, -1, dom.dr1)
                # println("West, 2nd order, c=$c")
            else
                # dϕdr = ( ϕ[c] - ϕ[c+wci])*dom.dr1
                # dv_r = -vf[c+wci, :] .*dom.dr1
                # pr =  dom.dr1

                # @. num += dϕdr*dv_r
                # den += dϕdr * pr
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, wci, -1, dom.dr1)
                # println("West, 1st order, c=$c")
            end
        # else
        #     println("No r stencil, c=$c")
        end
        
        # Z direction
        if usecell(acc, c+nci) # Use east
            if usecell(acc, c+2nci) # Second order north
                # dϕdz = (-1.5ϕ[c] + 2ϕ[c+nci] - 0.5ϕ[c+n⁺ci])*dom.dz1
                # pz = -1.5dom.dz1 
                # dv_z = (2vf[c+nci, :] .-0.5vf[c+n⁺ci,:]) .*dom.dz1

                # @. num += dϕdz*dv_z
                # den += dϕdz * pz
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, nci, 1, dom.dz1)
            else
                # dϕdz = (-ϕ[c] + ϕ[c+nci])*dom.dz1
                # pz = -dom.dz1
                # dv_z = vf[c+nci, :] .*dom.dz1

                # @. num += dϕdz*dv_z
                # den += dϕdz * pz
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, nci, 1, dom.dz1)
            end
        elseif usecell(acc, c+sci)
            if usecell(acc, c+2sci) # Second order south
                # dϕdz = ( 1.5ϕ[c] - 2ϕ[c+sci] + 0.5ϕ[c+s⁻ci])*dom.dz1
                # pz =     1.5dom.dz1 
                # dv_z = (          -2vf[c+sci, :] .+0.5vf[c+s⁻ci,:]).*dom.dz1

                # @. num += dϕdz*dv_z
                # den += dϕdz * pz
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, sci, -1, dom.dz1)
            else
                # dϕdz = ( ϕ[c] - ϕ[c+sci])*dom.dz1
                # pz =  dom.dz1
                # dv_z = -vf[c+sci, :] .*dom.dz1

                # @. num += dϕdz*dv_z
                # den += dϕdz * pz
                # println("c=$c")
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, sci, -1, dom.dz1)
            end
        end

        if den[1] == 0
            # pϕ= ϕ[c]
            # eϕ = checkbounds(Bool, acc, c+eci) ? ϕ[c+eci] : nothing
            # nϕ = checkbounds(Bool, acc, c+nci) ? ϕ[c+nci] : nothing
            # wϕ = checkbounds(Bool, acc, c+wci) ? ϕ[c+wci] : nothing
            # sϕ = checkbounds(Bool, acc, c+sci) ? ϕ[c+sci] : nothing
            # @warn "No identified stencil" c  pϕ eϕ nϕ wϕ sϕ
            # @debug "No identified stencil in fastmarch" c  
        else
            @. vf[c, :] = -num / den
        end

        # Add cell to accepted list
        acc[c] = true
        # println("Accepted: $(sum(acc))")

    end
    # debug
end

"""
    extrap_v_fastmarch(ϕ, T, dom::Domain, T_params)

Compute an extrapolated velocity field from T and ϕ.

Internally calls `compute_frontvel_withT` on positive half of Γ.
Using fast marching, instead of the PDE-based approach, to get second order accuracy more easily.

TODO: improve performance. Currently makes a lot of allocations, I think.
"""
function extrap_v_fastmarch(ϕ, T, dom::Domain, params)
    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Bf = identify_B(Γ, dom)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
    ϕ⁻ = ϕ .<= 0
    B⁻ = findall(Bf .& ϕ⁻)
    B⁺ = findall((ϕ⁻ .⊽ Γf) .& Bf) # Exclude Γ⁺

    vf = fill(0.0, dom.nr, dom.nz, 2)

    Qice = compute_Qice(ϕ, dom, params)
    icesurf = compute_icesurf(ϕ, dom)
    Qice_per_surf = Qice / icesurf

    # Accepted set
    acc = fill(false, dom.nr, dom.nz)


    # First, compute velocity on Γ⁺
    for c in Γ⁺
        vf[c, :] .= compute_frontvel_withT(T, ϕ, Tuple(c)..., dom, params, Qice_per_surf)
        acc[c] = true
    end

    # Second, fastmarch in positive half of B
    if isnan(sum(vf))
        println("Before fastmarch")
    end

    fastmarch_v!(vf, acc, B⁺, ϕ, dom)
    if isnan(sum(vf))
        println("Middle of fastmarch")
    end

    # Finally, fastmarch in negative half of B
    fastmarch_v!(vf, acc, B⁻, ϕ, dom)
    if isnan(sum(vf))
        println("After fastmarch")
    end

    vf
end


"""
"""
function extrap_v_pde(ϕi, Ti, dom, params)
    prop_t = 1.0
    # Precompute for velocity
    Qice = compute_Qice(ϕi, dom, params)
    icesurf = compute_icesurf(ϕi, dom)
    Qice_surf = Qice / icesurf
    
    Bf = identify_B(ϕi, dom)
    frontfunc(ir, iz) = compute_frontvel_withT(Ti, ϕi, ir, iz, dom, params, Qice_surf)
    vf = vector_extrap_from_front(ϕi, Bf, frontfunc, dom, prop_t)
    return vf
end