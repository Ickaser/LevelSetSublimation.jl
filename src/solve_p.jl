# export solve_T, solve_T_original
export solve_p_given_b
export solve_p, eval_b

"""
    solve_p(u, T, dom::Domain, params[, p0]; maxit=20, reltol=1e-6) where G<:AbstractArray

Iteratively compute the pressure profile for system state `u` and `T`.

There is a weak nonlinearity in the system, since the mass transfer coefficient `b` depends partially on pressure.
To treat this, use a guessed pressure `p0` (which, if not provided, is set everywhere to chamber pressure) to compute `b`,
then perform a linear solve for `p` using `solve_p_given_b`. At this point, recompute `b`, then recompute `p`.

In preliminary testing, this usually converges within 5 or 10 iterations.
Usually if it doesn't converge, it is because temperatures are outside the expected range, yielding crazy sublimation pressures.
(Occasionally it means I incorrectly wrote a finite difference somewhere.)

"""
function solve_p(u, Tf, T, dom::Domain, params; kwargs...) 
    # ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    ϕ = ϕ_T_from_u(u, dom)[1]
    # b = sum(eval_b(meanT, 0, dom, params))/dom.ntot
    b = eval_b(T, params[:p_ch], params)
    p0 = solve_p_given_b(ϕ, b, Tf, dom, params)
    solve_p(u, Tf, T, dom::Domain, params, p0; kwargs...)
end
function solve_p(u, Tf, T, dom::Domain, params, p0; maxit=20, reltol=1e-6) 
    # ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    ϕ = ϕ_T_from_u(u, dom)[1]

    relerr::Float64 = 0.0
    p⁺ = copy(p0)
    # Iterate up to maxit times
    for i in 1:maxit
        b = eval_b(T, p0, params)
        p⁺ .= solve_p_given_b(ϕ, b, Tf, dom, params)
        relerr = norm((p⁺ .- p0) ./ p⁺, Inf)
        if relerr < reltol
            # @info "Number of p iterations: $i"
            return p⁺
        end
        p0 .= p⁺
    end
    @info "Reached maximum iterations in p:" relerr maxit p⁺ T 
    return p⁺
end

"""
    eval_b(T, p, params)

Compute transport coefficient `b` as a function of space, given `T`, `p`, and `params`.

Uses the following fields from `params`, as documented in `sim_from_dict`: `\eps
`params` should have the following fields: `l`, `κ`, `R`, `Mw`, `μ`. 
`l`, and `κ` may be passed as scalars (and assumed as spatially uniform) or arrays (describing value throughout space, should match `Domain` dimensions).

When the simulation is started, all these values are converted to SI units and passed accordingly, 
so in practice there are no units to track.

If `κ=0`, no spatial variation due to pressure occurs.
"""
function eval_b(T, p, params)
    @unpack l, κ, R, Mw, μ = params
    b = @. Mw/R/T * (l*sqrt(R*T/Mw) + κ/μ*p)
end

"""
    eval_b_loc(T, p, ir, iz, params)

Locally evaluate transport coefficient (indexes into spatially varying `l` and `κ` if necessary).
"""
function eval_b_loc(T, p, ir, iz, params)
    @unpack l, κ, R, Mw, μ = params
    lloc = (length(l) > 1) ? l[ir,iz] : l
    κloc = (length(κ) > 1) ? κ[ir,iz] : κ
    b = Mw/R/T[ir,iz] * (lloc*sqrt(R*T[ir,iz]/Mw) + κloc/μ*p[ir,iz])
end

"""
    solve_p_given_b(ϕ, b, Tf, dom::Domain, params)

Compute 2D axisymmetric pseudosteady pressure profile for given values of level set function `ϕ`, temperature `T`, and transport coefficient `b`.

`b` is a dusty-gas transport coefficient for the pressure, which can vary spatially.
Homogeneous Neumann boundary conditions at `r=0`, `r=R`, `z=0`; Dirichlet on zero-level set (`p=p_sub`), Robin at top (`dp/dz = (p_ch-p)/Rp0`).
`params` should have fields: 
- `Rp0` : zero-thickness resistance offset, often written R0 in lyo literature
- `p_ch` : chamber (or vial) pressure at top surface

This implementation uses second-order finite differences, with linear extrapolation into Ω⁻.  
(For details, see [gibouFourthOrderAccurate2005](@cite).)  
Coefficients are all hard-coded here, unfortunately.
Neumann boundaries use a ghost point & BC to define ghost cell, then use same stencil as normal.
Coefficients computed in `gfm_extrap.ipynb`, using Sympy.  
"""
function solve_p_given_b(ϕ, b, Tf, dom::Domain, params) 
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Rp0, p_ch = params

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
            rhs[imx] = calc_psub(Tf[ir])
            continue
        end

        # Get local r values for use in Laplacian
        r = rgrid[ir]
        r1 = 1/r
        
        # θ_thresh = dr/dom.rmax
        θ_thresh = 0.05

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
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                psub_l = calc_psub(Tf_loc)
                if θr >= θ_thresh
                    pc += -2bp*dr2/θr
                    rhs[imx] -= 2bp*psub_l*dr2/θr #+ BC1*(r1 - 2dr1)
                else # Front is within dr/r cells of boundary
                    # Have an exact value given by BC + front
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                    # rhs[imx] = p_sub - θr*BC1*dr
                    # continue
                    pc += -2bp*dr2/(θr+1)
                    rhs[imx] -= 2psub_l*bp*dr2/(θr+1) 
                end
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
                # p. 119 of project notes
                θr = pϕ/(pϕ-wϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                psub_l = calc_psub(Tf_loc)
                if θr >= θ_thresh
                    pc += -2bp*dr2/θr
                    rhs[imx] -= 2psub_l*bp*dr2/θr + BC2*(r1 + 2dr1)
                else 
                    # First: use Neumann BC to get ghost cell left
                    # Second: extrapolate using left ghost cell across front
                    
                    pc += -2bp*dr2/(θr+1)
                    rhs[imx] -= 2psub_l*bp*dr2/(θr+1) + BC2*(r1 + 2dr1)
                end
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
            if eϕ <= 0 && wϕ <= 0
                # Pretend is in bulk, rather than treating with two ghost cells
                ec +=  bp*(1.0dr2 + 0.5dr1*r1) + dbr*0.5dr1
                pc += -bp*(2.0dr2)
                wc +=  bp*(1.0dr2 - 0.5dr1*r1) - dbr*0.5dr1
                rhs[imx] += 0

            elseif eϕ <= 0 # East ghost cell, across front
                θr = pϕ / (pϕ - eϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                psub_l = calc_psub(Tf_loc)
                if θr >= θ_thresh
                    # @info "hmm, east" θr dr ir iz
                    pc += (bp*(-(θr+1)*dr2 + (θr-1)*0.5dr1*r1) + dbr*(θr-1)*0.5dr1)/θr
                    wc +=  bp*(-0.5dr1*r1 + dr2) - dbr*0.5dr1 # Regular + gradient in b
                    rhs[imx] -= psub_l*(bp*(dr2+0.5dr1*r1) + dbr*0.5dr1)/θr # Dirichlet BC in ghost cell extrap
                else
                    # @info "hmm, east" θr dr ir iz
                    pc += -2bp*dr2 # Regular
                    wc += (bp*(2θr*dr2 - dr1*r1) - dbr*dr1)/(θr+1)
                    rhs[imx] -= psub_l*(bp*(2dr2+dr1*r1) + dbr*dr1)/(θr+1) # Dirichlet BC in ghost cell extrap
                end
            elseif wϕ <= 0 # West ghost cell across front
                θr = pϕ / (pϕ - wϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                psub_l = calc_psub(Tf_loc)
                if θr >= θ_thresh # Regular magnitude θ
                    pc += (bp*(-(θr+1)*dr2 + (1-θr)*0.5dr1*r1) + dbr*(1-θr)*0.5dr1)/θr
                    ec += bp*( 0.5dr1*r1 + dr2) + dbr*0.5dr1# Regular + b gradient
                    rhs[imx] -= psub_l*(bp*(dr2 - 0.5dr1*r1) - dbr*0.5dr1)/θr # Dirichlet BC in ghost cell extrap
                else # Very small θ
                    pc += -2bp*dr2 # Regular
                    ec += (bp*(2θr*dr2 + dr1*r1) + dbr*dr1)/(θr+1)
                    rhs[imx] -= psub_l*(bp*(2dr2 - dr1*r1) - dbr*dr1)/(θr+1) # Dirichlet BC in ghost cell extrap
                end

            else # Bulk, not at front 
                ec +=  bp*(1.0dr2 + 0.5dr1*r1) + dbr*0.5dr1
                pc += -bp*(2.0dr2)
                wc +=  bp*(1.0dr2 - 0.5dr1*r1) - dbr*0.5dr1
                rhs[imx] += 0
            end
        end

        # z direction discretization

        # For all z, psub is at given r
        psub_l = calc_psub(Tf[ir])
        # z direction boundaries
        if iz == 1
            # Zero flux BC
            BC3 = 0
            nϕ = ϕ[ir, iz+1]
            # Check for Stefan front
            if nϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θz = pϕ/(pϕ-nϕ)
                if θz > θ_thresh
                    pc += -2bp*dz2/θz
                    rhs[imx] -= 2*psub_l*bp*dz2/θz + 2*BC3*dz1
                else
                    # First: use Neumann BC to get ghost cell left
                    # Second: extrapolate using left ghost cell across front
                    pc += -2bp*dz2/(θz+1)
                    rhs[imx] -= 2psub_l*dz2*bp/(θz+1)
                end
            else
                # Using Neumann boundary to define ghost point: south T= north T - 2BC1*dr
                # p. 65, 66 of project notes
                # For gradient of conductivity: multiplied by gradient, which is BC, which is 0
                pc += -2bp*dz2
                nc +=  2bp*dz2
                rhs[imx] += -2*bp*BC3*dz1
            end
        elseif iz == nz
            # Robin BC

            sϕ = ϕ[ir, iz-1]
            dbz = (b[ir, iz]-b[ir, iz-1])*dz1
            # Check for Stefan front
            if sϕ < 0 # Front is within a cell of boundary
                # stefan_debug = true
                # p. 65 of project notes
                θz = pϕ/(pϕ-sϕ)
                if θz > θ_thresh
                    pc += -2bp*dz2/θz -2/Rp0*dz1 - dbz/bp/Rp0
                    rhs[imx] -= 2*psub_l*bp*dz2/θz + p_ch/Rp0*(dbz/bp + 2dz1)
                else
                    # First: use Robin BC to get ghost cell left
                    # Second: extrapolate using left ghost cell across front
                    pc += (-2bp*dz2 + dbz*dz1 - (dbz/bp + 2*θz*dz1)/Rp0)/(θz+1)
                    rhs[imx] -= (psub_l*(2dz2*bp - dbz*dz1) + p_ch/Rp0*(dbz/bp + 2θz*dz1))/(θz+1)
                end
                # No way to treat other equations, so cut it here
            else
                # Using Robin boundary to define ghost point: north p= south T - 2dz/bp*Rp0*(pi - p_ch)
                # p. 65, 66 of project notes
                pc += -2bp*dz2 - (2*dz1 + dbz/bp)/Rp0
                sc +=  2bp*dz2
                rhs[imx] -= p_ch*(dbz/bp + 2dz1)/Rp0
            end

        else # Bulk in z, still need to check for Stefan front
            nϕ = ϕ[ir, iz+1]
            sϕ = ϕ[ir, iz-1]
            dbz = (b[ir,iz+1]-b[ir, iz-1])*0.5*dz1
            if nϕ <= 0 && sϕ <= 0
                # Pretend is in bulk, rather than two ghost cells
                sc +=  bp*(1.0dz2) - dbz*0.5dz1
                pc += -bp*(2.0dz2) 
                nc +=  bp*(1.0dz2) + dbz*0.5dz1
                rhs[imx] += 0
            elseif nϕ <= 0 # North ghost cell
                θz = pϕ / (pϕ - nϕ)
                if θz >= θ_thresh
                    pc += (dbz*(θz-1)*dz1 - bp*(θz+1)*dz2 )/θz
                    sc += bp*dz2 - dbz*0.5*dz1
                    rhs[imx] -= psub_l*(bp*dz2 + dbz*0.5*dz1)/θz
                else
                    pc += -2bp*dz2
                    sc += (2*bp*θz*dz2 - dbz*dz1)/(θz+1)
                    rhs[imx] -= psub_l*(2bp*dz2 + dbz*dz1)/(θz+1)
                end
            elseif sϕ <= 0
                θz = pϕ / (pϕ - sϕ)
                if θz >= θ_thresh
                    pc += (-bp*dz2*(θz+1) + dbz*(1-θz)*dz1)/θz
                    nc += bp*dz2 + dbz*0.5dz1
                    rhs[imx] -= psub_l*(bp*dz2 - dbz*0.5dz1)/θz
                else
                    pc += -2bp*dz2
                    nc += (2*bp*θz*dz2 + dbz*dz1)/(θz+1)
                    rhs[imx] -= psub_l*(2bp*dz2 - dbz*dz1)/(θz+1)
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
