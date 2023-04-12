export fastmarch_v!, extrap_v_fastmarch

export compute_frontvel_1



    
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
    extrap_v_fastmarch(ϕ, T, p, dom::Domain, params)

Compute an extrapolated velocity field from `ϕ`, `T`, and `p`.

Internally calls `compute_frontvel_withT` on positive half of Γ.
Using fast marching, instead of the PDE-based approach, to get second order accuracy more easily.

TODO: improve performance. Currently makes a lot of allocations, I think.
"""
function extrap_v_fastmarch(u, T, p, dom::Domain, params)
    ϕ = ϕ_T_from_u(u, dom)[1]
    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
    ϕ⁻ = ϕ .<= 0
    # Bf = identify_B(Γ, dom)
    # B⁻ = findall(Bf .& ϕ⁻)
    # B⁺ = findall((ϕ⁻ .⊽ Γf) .& Bf) # Exclude Γ⁺
    Ω⁻ = findall(ϕ⁻)
    Ω⁺ = findall(ϕ⁻ .⊽ Γf ) # Exclude Γ⁺

    # Accepted set
    acc = fill(false, dom.nr, dom.nz)

    # First, compute velocity on Γ⁺



    vf, dϕdx_all = compute_frontvel_mass(u, T, p, dom, params)

    # pl1 = heat(vf[:,:,1], dom)
    # plot_contour(ϕ, dom)
    acc[Γ⁺] .= true

    # Second, fastmarch in positive half of Ω
    fastmarch_v!(vf, acc, Ω⁺, ϕ, dom)
    # pl2 = heat(vf[:,:,1], dom)

    # Finally, fastmarch in negative half of Ω
    fastmarch_v!(vf, acc, Ω⁻, ϕ, dom)
    # pl3 = heat(vf[:,:,1], dom)

    # pl4 = heat(ϕ, dom)
    # plot_contour(ϕ, dom)

    # display(plot(pl1, pl2, pl3, pl4))
    # @info "extrap" argmin(vf[:,:,1])

    vf, dϕdx_all
end
