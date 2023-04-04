# export solve_T, solve_T_original
export solve_p_given_b
export solve_p

" Takes T as Kelvin, returns P in Pa"
function calc_psub(T)
    ai = [-0.212144006e2,  0.273203819e2,  -0.610598130e1]
    bi = [0.333333333e-2,  0.120666667e1,  0.170333333e1]
    θ = T/275.16
    lnπ = sum(ai .* θ .^bi) / θ
    exp(lnπ)*611.657
end

"""
    eval_b(T, p, params)

Compute transport coefficient `b` as a function of space, given `T`, `p`, and `params`.

`params` should have the following fields:
- `ϵ` : porosity of porous medium
- `l` : dusty gas model constant: characteristic length for Knudsen diffusion
- `κ` : dusty gas model constant: length^2 corresponding loosely to Darcy's Law permeability
- `R` : universal gas constant, with appropriate units 
- `Mw`: molecular weight of species (water), with appropriate units
- `μ` : dynamic viscosity of species (water), with appropriate units
`ϵ`, `l`, and `κ` may be passed as scalars (and assumed as spatially uniform) or arrays (describing value throughout space, should match `Domain` dimensions).

If `κ=0`, no spatial variation due to pressure occurs.
"""
function eval_b(T, p, params)
    @unpack ϵ, l, κ, R, Mw, μ = params
    b = @. 1/R/T/Mw * (l*sqrt(R*T/Mw) + κ/μ*p)
end

"""
    solve_p_given_b(ϕ, b, dom::Domain, params)

Compute 2D axisymmetric pseudosteady pressure profile for given values of level set function `ϕ`, temperature `T`, and transport coefficient `b`.

`b` is a dusty-gas transport coefficient for the pressure. 
- If `b` is a scalar, will be treated as spatially constant (which assumes `κ = 0`). 
- If `b` an array of same size as `ϕ` and `T`, gradients will be computed.
Homogeneous Neumann boundary conditions at `r=0`, `r=R`, `z=0`; Dirichlet on zero-level set (`p=p_sub`) and top (`p=p_ch`).
`params` should have fields: 
- `p_ch` : chamber (or vial) pressure at top surface
- `p_sub` : sublimation pressure on interface

This implementation uses second-order finite differences, with linear extrapolation into Ω⁻.  
(For details, see [gibouFourthOrderAccurate2005](@cite).)  
Coefficients are all hard-coded here, unfortunately.
Neumann boundaries use a ghost point & BC to define ghost cell, then use same stencil as normal.
Coefficients computed in `gfm_extrap.ipynb`, using Sympy.  
"""
# First: version with no spatial variation
function solve_p_given_b(ϕ, b::Float64, dom::Domain, params)
    @info "Using scalar b for p"
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack p_ch, p_sub = params

    # To prevent blowup, artificially add some corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    if minimum(ϕ) > 0
        # @info "Solving heat equation without any ice, artificially introducing some"
        ϕ = copy(ϕ)
        ϕ[argmin(ϕ)] = - max(dr, dz)
    end

    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = Vector{Float64}(undef, 0)
    vcr = (vals, cols, rows)
    rhs = fill(0.0, ntot)

    for iz in 1:nz, ir in 1:nr
        # Row position in matrix: r is small iteration, z is outer iteration
        imx = ir + (iz-1)*nr

        pϕ = ϕ[ir, iz]

        # Check if in frozen domain; if so, fix pressure at 1.1psub
        if pϕ <= 0
            add_to_vcr!(vcr, dom, imx, (0, 0), 1)
            # rhs[imx] = 1.1p_sub
            rhs[imx] = p_sub
            continue
        end

        # Get local r values for use in Laplacian
        r = rgrid[ir]
        r1 = 1/r

        # Stencil values: initialize to 0
        ec = 0
        pc = 0
        wc = 0
        sc = 0
        nc = 0

        # R direction boundaries
        if ir == 1
            # Symmetry BC
            BC1 = 0
            eϕ = ϕ[ir+1, iz]
            # Check for Stefan front
            if eϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θr = pϕ/(pϕ-eϕ)
                # Have an exact value given by BC + front
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub - θr*BC1*dr
                continue
                # No way to treat other equations, so cut it here
            else
                # Using Neumann boundary to define ghost point: west T= east T - 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2dr2
                ec +=  2dr2
                # rhs[imx] += BC1*(2*dr1 - r1)
                rhs[imx] += 0 # r=0, so 1/r = NaN
            end
        elseif ir == nr
            # Zero flux BC
            BC2 = 0
            wϕ = ϕ[ir-1, iz]
            if wϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θr = pϕ/(pϕ-wϕ)
                # Have an exact value given by BC + front
                # No way to treat other equations, so add to matrix and stop here
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub + θr*BC2*dr
                continue
            else
                # Using Neumann boundary to define ghost point: east T= west T + 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2dr2
                wc +=  2dr2
                rhs[imx] += BC2*(-2*dr1 - r1)
            end

        else # r direction bulk
            # Check for Stefan boundary
            eϕ = ϕ[ir+1, iz]
            wϕ = ϕ[ir-1, iz]
            if eϕ <= 0 # East ghost cell, across front
                θr = pϕ / (pϕ - eϕ)
                if θr >= dr
                    pc += -2b*dr2 # Regular 
                    pc += b*(θr-1)/θr*(0.5dr1*r1+ dr2) # Due to ghost cell extrapolation
                    wc += b*(-0.5dr1*r1 + dr2) # Regular
                    rhs[imx] -= p_sub*b*(0.5*dr+r) *dr2 *r1/θr # Dirichlet BC in ghost cell extrap
                else
                    pc += -2b*dr2 # Regular
                    wc += b*(-0.5dr1*r1 + dr2) # Regular
                    wc += b*(dr+2r)*(θr-1)*0.5dr2*r1/(θr+1) # Due to ghost cell extrapolation 
                end
            elseif wϕ <= 0 # West ghost cell across front
                θr = pϕ / (pϕ - wϕ)
                if θr >= dr # Regular magnitude θ
                    pc += -2b*dr2 # Regular 
                    pc += b*(-dr+2r)*(θr-1)*0.5dr2*r1/θr # Due to ghost cell extrapolation
                    ec += b*( 0.5dr1*r1 + dr2) # Regular
                    rhs[imx] -= p_sub*b*(-0.5dr+r) *dr2 *r1/θr # Dirichlet BC in ghost cell extrap
                else # Very small θ
                    pc += -2b*dr2 # Regular
                    ec += b*( 0.5dr1*r1 + dr2) # Regular
                    ec += b*(-dr+2r)*(θr-1)/(θr+1)*0.5dr2*r1 # Due to ghost cell extrapolation
                    rhs[imx] -= p_sub*b*(-dr+2r) *dr2 *r1/(θr+1) # Dirichlet BC in ghost cell extrap
                end

            else # Bulk, not at front 
                ec +=  1.0dr2 + 0.5dr1*r1
                pc += -2.0dr2
                wc +=  1.0dr2 - 0.5dr1*r1
                rhs[imx] += 0
            end
        end

        # z direction discretization

        # z direction boundaries
        if iz == 1
            # Zero flux BC
            BC3 = 0
            nϕ = ϕ[ir, iz+1]
            # Check for Stefan front
            if nϕ < 0 # Front is within a cell of boundary
                # stefan_debug = true
                # p. 65 of project notes
                θz = pϕ/(pϕ-nϕ)
                # Have an exact value given by BC + front
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub + θz*BC3*dz
                continue
                # No way to treat other equations, so cut it here
            else
                # Using Neumann boundary to define ghost point: south T= north T - 2BC1*dr
                # p. 65, 66 of project notes
                pc += -2dz2
                nc +=  2dz2
                rhs[imx] += -2*BC3*dz1
            end
        elseif iz == nz
            # Dirichlet BC
            BC4 = p_ch
            add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
            rhs[imx] = BC4
            continue

        else # Bulk in z, still need to check for Stefan front
            nϕ = ϕ[ir, iz+1]
            sϕ = ϕ[ir, iz-1]
            if nϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - nϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, north")
                if θz > dz
                    pc += -b*dz2*(θz+1)/θz
                    sc += b*dz2
                    rhs[imx] -= p_sub*b*dz2/θz
                else
                    pc += -2b*dz2
                    sc += 2*b*θz*dz2/(θz+1)
                    rhs[imx] -= 2p_sub*b*dz2/(θz+1)
                end
            elseif sϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - sϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, south")
                if θz > dz
                    pc += -b*dz2*(θz+1)/θz
                    nc += b*dz2
                    rhs[imx] -= p_sub*b*dz2/θz
                else
                    pc += -2b*dz2
                    nc += 2*b*θz*dz2/(θz+1)
                    rhs[imx] -= 2p_sub*b*dz2/(θz+1)
                end

            else # Bulk, no Stefan front
                sc +=  1.0dz2
                pc += -2.0dz2
                nc +=  1.0dz2
                rhs[imx] += 0
            end
        end

        # Assign all computed stencil values into matrix
        pc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
        ec != 0 && add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
        wc != 0 && add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
        nc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
        sc != 0 && add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)

    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    sol = mat_lhs \ rhs
    psol = reshape(sol, nr, nz)
end

# Version with spatial variation
function solve_p_given_b(ϕ, b::T, dom::Domain, params) where T<:AbstractArray
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack p_ch, p_sub = params

    # To prevent blowup, artificially add some corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    if minimum(ϕ) > 0
        # @info "Solving heat equation without any ice, artificially introducing some"
        ϕ = copy(ϕ)
        ϕ[argmin(ϕ)] = - max(dr, dz)
    end

    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = Vector{Float64}(undef, 0)
    vcr = (vals, cols, rows)
    rhs = fill(0.0, ntot)

    for iz in 1:nz, ir in 1:nr
        # Row position in matrix: r is small iteration, z is outer iteration
        imx = ir + (iz-1)*nr

        pϕ = ϕ[ir, iz]
        bp = b[ir, iz]

        # Check if in frozen domain; if so, fix pressure at arbitrary value
        if pϕ <= 0
            add_to_vcr!(vcr, dom, imx, (0, 0), 1)
            # rhs[imx] = 1.1p_sub
            rhs[imx] = p_sub
            # rhs[imx] = p_ch
            continue
        end

        # Get local r values for use in Laplacian
        r = rgrid[ir]
        r1 = 1/r

        # Stencil values: initialize to 0
        ec = 0
        pc = 0
        wc = 0
        sc = 0
        nc = 0

        # R direction boundaries
        if ir == 1
            # Symmetry BC
            BC1 = 0
            eϕ = ϕ[ir+1, iz]
            # Check for Stefan front
            if eϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θr = pϕ/(pϕ-eϕ)
                # Have an exact value given by BC + front
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub - θr*BC1*dr
                continue
                # No way to treat other equations, so cut it here
            else
                # Using Neumann boundary to define ghost point: west T= east T - 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                # For gradient of conductivity: multiplied by gradient, which is BC, which is 0
                pc += -2bp*dr2
                ec +=  2bp*dr2
                # rhs[imx] += BC1*(2*dr1 - r1)
                rhs[imx] += 0 # r=0, so 1/r = NaN
            end
        elseif ir == nr
            # Zero flux BC
            BC2 = 0
            wϕ = ϕ[ir-1, iz]
            if wϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θr = pϕ/(pϕ-wϕ)
                # Have an exact value given by BC + front
                # No way to treat other equations, so add to matrix and stop here
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub + θr*BC2*dr
                continue
            else
                # Using Neumann boundary to define ghost point: east T= west T + 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                # For gradient of conductivity: multiplied by gradient, which is BC, which is 0
                pc += -2bp*dr2 + 0 
                wc +=  2bp*dr2 + 0
                rhs[imx] += 0
            end

        else # r direction bulk
            # Check for Stefan boundary
            dbr = (b[ir+1, iz]-b[ir-1, iz])*0.5*dr1
            eϕ = ϕ[ir+1, iz]
            wϕ = ϕ[ir-1, iz]
            if eϕ <= 0 # East ghost cell, across front
                θr = pϕ / (pϕ - eϕ)
                if θr >= dr
                    pc += (bp*(-(θr+1)*dr2 + (θr-1)*0.5dr1*r1) + dbr*(θr-1)*0.5dr1)/(θr+1)
                    wc += bp*(0.5dr1*r1 + dr2) - dbr*0.5dr1 # Regular + gradient in b
                    rhs[imx] -= p_sub*(bp*(dr2+0.5dr1*r1) + dbr*0.5dr1)/θr # Dirichlet BC in ghost cell extrap
                else
                    # @info "hmm, east" θr dr ir iz
                    pc += -2bp*dr2 # Regular
                    wc += (bp*(2dr2 - dr1*r1) - dbr*dr1)/(θr+1)
                    rhs[imx] -= p_sub*(bp*(2dr2+dr1*r1) + dbr*dr1)/(θr+1) # Dirichlet BC in ghost cell extrap
                end
            elseif wϕ <= 0 # West ghost cell across front
                θr = pϕ / (pϕ - wϕ)
                if θr >= dr # Regular magnitude θ
                    pc += (bp*(-(θr+1)*dr2 + (1-θr)*0.5dr1*r1) + dbr*(1-θr)*0.5dr1)/θr
                    ec += bp*( 0.5dr1*r1 + dr2) + dbr*0.5dr1# Regular + b gradient
                    rhs[imx] -= p_sub*(bp*(dr2 - 0.5dr1*r1) - dbr*0.5dr1)/θr # Dirichlet BC in ghost cell extrap
                else # Very small θ
                    # @info "hmm, west" θr dr ir iz
                    pc += -2bp*dr2 # Regular
                    ec += (bp*(dr1*r1 + 2θr*dr2) + dbr*dr1)/(θr+1)
                    rhs[imx] -= p_sub*(bp*(2dr2-dr1*r1)-dbr*dr1)/(θr+1) # Dirichlet BC in ghost cell extrap
                end

            else # Bulk, not at front 
                ec +=  bp*(1.0dr2 + 0.5dr1*r1) + dbr*0.5dr1
                pc += -bp*(2.0dr2)
                wc +=  bp*(1.0dr2 - 0.5dr1*r1) - dbr*0.5dr1
                rhs[imx] += 0
            end
        end

        # z direction discretization

        # z direction boundaries
        if iz == 1
            # Zero flux BC
            BC3 = 0
            nϕ = ϕ[ir, iz+1]
            # Check for Stefan front
            if nϕ < 0 # Front is within a cell of boundary
                # stefan_debug = true
                # p. 65 of project notes
                θz = pϕ/(pϕ-nϕ)
                # Have an exact value given by BC + front
                add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                rhs[imx] = p_sub + θz*BC3*dz
                continue
                # No way to treat other equations, so cut it here
            else
                # Using Neumann boundary to define ghost point: south T= north T - 2BC1*dr
                # p. 65, 66 of project notes
                # For gradient of conductivity: multiplied by gradient, which is BC, which is 0
                pc += -2bp*dz2
                nc +=  2bp*dz2
                rhs[imx] += -2*bp*BC3*dz1
            end
        elseif iz == nz
            # Dirichlet BC
            BC4 = p_ch
            add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
            rhs[imx] = BC4
            continue

        else # Bulk in z, still need to check for Stefan front
            nϕ = ϕ[ir, iz+1]
            sϕ = ϕ[ir, iz-1]
            dbz = (b[iz+1]-b[iz-1])*0.5*dz1
            if nϕ <= 0 # North ghost cell
                θz = pϕ / (pϕ - nϕ)
                if θz > dz
                    pc += (0.5dbz*(θz-1)*dz1 - bp*(θz+1)*dz2 )/θz
                    sc += bp*dz2 - dbz*0.5*dz1
                    rhs[imx] -= p_sub*(bp*dz2 + dbz*0.5*dz1)/θz
                else
                    # @info "hmm, north" θz dz ir iz
                    pc += -2bp*dz2
                    sc += (2*bp*θz*dz2 - dbz*dz1)/(θz+1)
                    rhs[imx] -= p_sub*(2bp*dz2 + dbz*dz1)/(θz+1)
                end
            elseif sϕ <= 0
                θz = pϕ / (pϕ - sϕ)
                if θz > dz
                    pc += (-bp*dz2*(θz+1) + dbz*(1-θz)*dz1)/θz
                    nc += bp*dz2 + dbz*0.5dz1
                    rhs[imx] -= p_sub*(bp*dz2 - dbz*0.5dz1)/θz
                else
                    # @info "hmm, south" θz dz ir iz
                    pc += -2bp*dz2
                    nc += (2*bp*θz*dz2 + dbz*dz1)/(θz+1)
                    rhs[imx] -= p_sub*(2bp*dz2 - dbz*dz1)/(θz+1)
                end

            else # Bulk, no Stefan front
                sc +=  bp*(1.0dz2) - dbz*0.5dz1
                pc += -bp*(2.0dz2) 
                nc +=  bp*(1.0dz2) + dbz*0.5dz1
                rhs[imx] += 0
            end
        end

        # Assign all computed stencil values into matrix
        pc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
        ec != 0 && add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
        wc != 0 && add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
        nc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
        sc != 0 && add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)

    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    sol = mat_lhs \ rhs
    psol = reshape(sol, nr, nz)
end

function solve_p(ϕ, T, dom::Domain, params; p0::Union{Nothing, G}=nothing, maxit=10, reltol=1e-6) where G<:AbstractArray
    if p0 === nothing
        # meanT = sum(T[ϕ .>0]) / sum(ϕ .> 0)
        # b = sum(eval_b(meanT, 0, dom, params))/dom.ntot
        b = eval_b(T, 0, params)
        p0 = solve_p_given_b(ϕ, b, dom, params)
    end

    relerr::Float64 = 0.0
    p⁺ = copy(p0)
    # Iterate 5 times
    for i in 1:maxit
        b = eval_b(T, p0, params)
        p⁺ = solve_p_given_b(ϕ, b, dom, params)
        relerr = maximum(abs.(p⁺ .- p0) ./ p⁺)
        # @info "Maximum relative error in p after $i iterations:" relerr
        if relerr < reltol
            return p⁺
        end
        p0 = p⁺
    end
    @info "Reached maximum iterations in p:" relerr maxit
    return p⁺
end