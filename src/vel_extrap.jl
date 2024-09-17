export fastmarch_v!, extrap_v_fastmarch!

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
        # num = fill(0.0, 2)
        # den = [0.0]
        num = zeros(eltype(ϕ), 2)
        den = zeros(eltype(ϕ), 1)

        # R direction
        if usecell(acc, c+eci) # Use east
            if usecell(acc, c+2eci) # Second order east
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, eci, 1, dom.dr1)
            else
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, eci, 1, dom.dr1)
            end
        elseif usecell(acc, c+wci)
            if usecell(acc, c+2wci) # Second order west
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, wci, -1, dom.dr1)
            else
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, wci, -1, dom.dr1)
            end
        end
        
        # Z direction
        if usecell(acc, c+nci) # Use north
            if usecell(acc, c+2nci) # Second order north
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, nci, 1, dom.dz1)
            else
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, nci, 1, dom.dz1)
            end
        elseif usecell(acc, c+sci)
            if usecell(acc, c+2sci) # Second order south
                ddx_fastmarch_2nd!(num, den, ϕ, vf, c, sci, -1, dom.dz1)
            else
                ddx_fastmarch_1st!(num, den, ϕ, vf, c, sci, -1, dom.dz1)
            end
        end

        if den[1] == 0
            @debug "No identified stencil in fastmarch" c  
        else
            @. vf[c, :] = -num / den
        end
        # Add cell to accepted list
        acc[c] = true
    end
end

"""
    extrap_v_fastmarch(v_front, u, dom::Domain)

Compute an extrapolated velocity field in place. `v_front` should be nonzero only on Γ⁺, positive side of interface.

Using fast marching, instead of the PDE-based approach, to get second order accuracy more easily.
"""
function extrap_v_fastmarch!(v_front, u, dom::Domain)
    ϕ = reshape(u[iϕ(dom)], size(dom))
    Γf = identify_Γ(ϕ, dom)
    Γ = findall(Γf)
    Γ⁺ = [c for c in Γ if ϕ[c]>0]
    ϕ⁻ = ϕ .<= 0
    Ω⁻ = findall(ϕ⁻)
    Ω⁺ = findall(ϕ⁻ .⊽ Γf ) # Exclude Γ⁺, ⊽ is "nor" operator

    # Accepted set
    acc = fill(false, dom.nr, dom.nz)
    acc[Γ⁺] .= true
    if extrema(v_front[.~ acc, :]) != (0.0, 0.0) # Everywhere but "accepted" region, velocity is 0
        @warn "Velocity field may not be computed on correct cells before extrapolation." acc v_front
    end

    # First, fastmarch in positive half of Ω
    fastmarch_v!(v_front, acc, Ω⁺, ϕ, dom)

    # Second, fastmarch in negative half of Ω
    fastmarch_v!(v_front, acc, Ω⁻, ϕ, dom)
    return nothing
end
