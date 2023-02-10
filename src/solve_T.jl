export solve_T, solve_T_original
# y


# Construct sparse array with rows, cols, vals format, then construct sparse matrix
"""
    solve_T_original(ϕ, dom::Domain, params)

Compute 2D axisymmetric temperature profile for given level set function `ϕ`.

`params` should have fields: `Q_gl`, `Q_sh`, `Q_ck`, `k`, `Tf`
This implementation uses second-order finite differences, with constant extrapolation into Ω⁻.
"""
function solve_T_original(ϕ, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Q_gl, Q_sh, Q_ck, k, Tf = params
    # To prevent blowup, artificially add some  corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    if minimum(ϕ) > 0
        @warn "Solving heat equation without any ice, artificially introducing some"
        ϕ = copy(ϕ)
        ϕ[1,2] = -max(dr, dz)
        ϕ[2,1] = -max(dr, dz)
    end
    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = Vector{Float64}(undef, 0)
    rhs = fill(0.0, ntot)

    for iz in 1:nz
        for ir in 1:nr
            # Get local r, z values
            r = rgrid[ir]
            r1 = 1/r
            z = zgrid[iz]
            
            # Row position in matrix: r is small iteration, z is outer iteration
            imx = ir + (iz-1)*nr

            # Check if in frozen domain; if so, fix temperature
            if ϕ[ir, iz] <= 0
                push!(vals, 1.0); push!(cols, imx); push!(rows, imx)
                rhs[imx] = Tf
                continue
            end
            #

            # Check if on Stefan boundary
            # (actually I shouldn't need to do this but this is where I would do it if I did)
            
            # Check if on boundaries or in bulk, set up BC or diffeq  

            # R direction boundaries
            if ir == 1
                # Symmetry BC
                push!(vals, -1.5dr1); push!(cols, imx  ); push!(rows, imx) # P cell
                push!(vals,  2.0dr1); push!(cols, imx+1); push!(rows, imx) # E cell
                push!(vals, -0.5dr1); push!(cols, imx+2); push!(rows, imx) # E+ cell

                rhs[imx] = 0.0
                # No z component, to enforce BC
                continue
            elseif ir == nr
                # Constant flux BC
                push!(vals,  0.5dr1); push!(cols, imx-2); push!(rows, imx) # W- cell
                push!(vals, -2.0dr1); push!(cols, imx-1); push!(rows, imx) # W cell
                push!(vals,  1.5dr1); push!(cols, imx  ); push!(rows, imx) # P cell

                rhs[imx] = Q_gl / k
                # No z component, to enforce BC
                continue
            end
            # z direction boundaries
            if iz == 1
                # # Constant flux BC
                push!(vals, -1.5dz1); push!(cols, imx  ); push!(rows, imx) # P cell
                push!(vals,  2.0dz1); push!(cols, imx+ nr); push!(rows, imx) # N cell
                push!(vals, -0.5dz1); push!(cols, imx+2nr); push!(rows, imx) # N+ cell

                rhs[imx] = -Q_sh / k

                # No r component, to enforce BC
                continue
            elseif iz == nz
                # Adiabatic BC
                push!(vals,  0.5dz1); push!(cols, imx-2nr); push!(rows, imx) # S- cell
                push!(vals, -2.0dz1); push!(cols, imx- nr); push!(rows, imx) # S cell
                push!(vals,  1.5dz1); push!(cols, imx  ); push!(rows, imx) # P cell

                rhs[imx] = 0.0
                # No r component, to enforce BC
                continue
            end


            # Bulk Eq: Cylindrical laplacian, second order finite diffs, d2/dr2 + 1/r d/dr + d2/dz2
            if ir >= 2 && ir <= nr-1 && iz >= 2 && iz <= nr-1
                push!(vals,                     1.0dz2); push!(cols, imx-nr); push!(rows, imx) # S cell
                push!(vals,  1.0dr2 - 0.5dr1*r1       ); push!(cols, imx- 1); push!(rows, imx) # W cell
                push!(vals, -2.0dr2            -2.0dz2); push!(cols, imx   ); push!(rows, imx) # P cell
                push!(vals,  1.0dr2 + 0.5dr1*r1       ); push!(cols, imx+ 1); push!(rows, imx) # E cell
                push!(vals,                     1.0dz2); push!(cols, imx+nr); push!(rows, imx) # N cell
                rhs[imx] = - Q_ck / k
            end
        end
    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    sol = mat_lhs \ rhs
    T = reshape(sol, nr, nz)
end
        

"""
    solve_T(ϕ, dom::Domain, params)

Compute 2D axisymmetric T profile, returning ghost cell values, for given level set function `ϕ`.

Neumann boundary conditions on all rectangular boundaries; Dirichlet on zero-level set.
`params` should have fields: 
- `Q_gl` 
- `Q_sh` 
- `Q_ck`
- `k` 
- `Tf`  
This implementation uses second-order finite differences, with linear extrapolation into Ω⁻.  
Coefficients are all hard-coded here, unfortunately.
(For more on extrapolation, see Gibou et al., 2002, "Second-Order-Accurate ... Poisson ... ")  
Coefficients computed in `gfm_extrap.ipynb`, using Sympy.  
(For higher order, see Gibou and Fedkiw, 2005, "A fourth order accurate discretization ... Laplace ... ")  
"""
function solve_T(ϕ, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Q_gl, Q_sh, Q_ck, k, Tf = params
    # To prevent blowup, artificially add some  corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    if minimum(ϕ) > 0
        @warn "Solving heat equation without any ice, artificially introducing some"
        ϕ = copy(ϕ)
        ϕ[1,2] = -max(dr, dz)
        ϕ[2,1] = -max(dr, dz)
    end
    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = Vector{Float64}(undef, 0)

    vcr = (vals, cols, rows)
    rhs = fill(0.0, ntot)

    for iz in 1:nz
        for ir in 1:nr
            # Get local r values for use in Laplacian
            r = rgrid[ir]
            r1 = 1/r
            # z = zgrid[iz]
            
            # Row position in matrix: r is small iteration, z is outer iteration
            imx = ir + (iz-1)*nr

            # Check if in frozen domain; if so, fix temperature
            if ϕ[ir, iz] <= 0
                # push!(vals, 1.0); push!(cols, imx); push!(rows, imx)
                add_to_vcr!(vcr, dom, imx, (0, 0), 1)
                rhs[imx] = Tf
                continue
            end
            #

            # Check if on boundaries or in bulk, set up BC or diffeq  
            # TODO 
            # Currently, boundaries are not checking for interface crossing.

            # R direction boundaries
            if ir == 1
                # Symmetry BC
                # push!(vals, -1.5dr1); push!(cols, imx  ); push!(rows, imx) # P cell
                # push!(vals,  2.0dr1); push!(cols, imx+1); push!(rows, imx) # E cell
                # push!(vals, -0.5dr1); push!(cols, imx+2); push!(rows, imx) # E+ cell
                add_to_vcr!(vcr, dom, imx, ( 0, 0), k* -1.5dr1) # P cell
                add_to_vcr!(vcr, dom, imx, ( 1, 0), k*  2.0dr1) # E cell
                add_to_vcr!(vcr, dom, imx, ( 2, 0), k* -0.5dr1) # E+ cell

                rhs[imx] = 0.0
                # No z component, to enforce BC
                continue
            elseif ir == nr
                # Constant flux BC
                # push!(vals,  0.5dr1); push!(cols, imx-2); push!(rows, imx) # W- cell
                # push!(vals, -2.0dr1); push!(cols, imx-1); push!(rows, imx) # W cell
                # push!(vals,  1.5dr1); push!(cols, imx  ); push!(rows, imx) # P cell
                add_to_vcr!(vcr, dom, imx, ( 0, 0), k*  1.5dr1) # P cell
                add_to_vcr!(vcr, dom, imx, (-1, 0), k* -2.0dr1) # W cell
                add_to_vcr!(vcr, dom, imx, (-2, 0), k*  0.5dr1) # W- cell

                rhs[imx] = Q_gl 
                # No z component, to enforce BC
                continue
            end
            # z direction boundaries
            if iz == 1
                # # Constant flux BC
                # push!(vals, -1.5dz1); push!(cols, imx  ); push!(rows, imx) # P cell
                # push!(vals,  2.0dz1); push!(cols, imx+ nr); push!(rows, imx) # N cell
                # push!(vals, -0.5dz1); push!(cols, imx+2nr); push!(rows, imx) # N+ cell
                add_to_vcr!(vcr, dom, imx, ( 0, 0), k* -1.5dz1) # P cell
                add_to_vcr!(vcr, dom, imx, ( 0, 1), k*  2.0dz1) # E cell
                add_to_vcr!(vcr, dom, imx, ( 0, 2), k* -0.5dz1) # E+ cell

                rhs[imx] = -Q_sh 

                # No r component, to enforce BC
                continue
            elseif iz == nz
                # Adiabatic BC
                # push!(vals,  0.5dz1); push!(cols, imx-2nr); push!(rows, imx) # S- cell
                # push!(vals, -2.0dz1); push!(cols, imx- nr); push!(rows, imx) # S cell
                # push!(vals,  1.5dz1); push!(cols, imx  ); push!(rows, imx) # P cell
                add_to_vcr!(vcr, dom, imx, ( 0, 0), k* 1.5dz1) # P cell
                add_to_vcr!(vcr, dom, imx, ( 0,-1), k*-2.0dz1) # E cell
                add_to_vcr!(vcr, dom, imx, ( 0,-2), k* 0.5dz1) # E+ cell

                rhs[imx] = 0.0
                # No r component, to enforce BC
                continue
            end


            ec = 0
            pc = 0
            wc = 0
            sc = 0
            nc = 0

            # Check if on Stefan boundary
            # if 
            # Conditions for Stefan boundary:
            stencil = [CartesianIndex(i) for i in [(0,0), (1, 0), (0, 1), (-1, 0), (0, -1)]] # p, e, n, w, s
            # println(stencil)
            self = CartesianIndex((ir, iz))
            loc_stencil = [s + self for s in stencil]
            # ϕ > 0 implies dried domain
            ϕst = ϕ[loc_stencil]
            e_front = (ϕst[2] < 0)
            n_front = (ϕst[3] < 0)
            w_front = (ϕst[4] < 0)
            s_front = (ϕst[5] < 0)
            r_bulk = !e_front && !w_front
            z_bulk = !n_front && !s_front

            if e_front
                θr = ϕst[1] / (ϕst[1] - ϕst[2])
                if θr >= dr
                    pc += -2k*dr2 # Regular 
                    pc += k*(θr-1)*(0.5dr1/r+ dr2)/θr # Due to ghost cell extrapolation
                    wc += k*(-0.5dr1/r + dr2) # Regular
                    rhs[imx] -= Tf*k*(0.5*dr+r) *dr2 /r/θr # Dirichlet BC in ghost cell extrap
                else
                    pc += -2k*dr2 # Regular
                    wc += k*(-0.5dr1/r + dr2) # Regular
                    wc += k*(dr+2r)*(θr-1)*0.5dr2/r/(θr+1) # Due to ghost cell extrapolation in extrapolation
                    rhs[imx] -= Tf*k*(0.5*dr+r) *dr2 /r/(θr+1) # Dirichlet BC in ghost cell extrap
                end
            elseif w_front
                θr = ϕst[1] / (ϕst[1] - ϕst[4])
                if θr >= dr
                    pc += -2k*dr2 # Regular 
                    pc += k*(-dr+2r)*(θr-1)*0.5dr2/r/θr # Due to ghost cell extrapolation
                    ec += k*( 0.5dr1/r + dr2) # Regular
                    rhs[imx] -= Tf*k*(-0.5dr+r) *dr2 /r/θr # Dirichlet BC in ghost cell extrap
                else
                    pc += -2k*dr2 # Regular
                    ec += k*( 0.5dr1/r + dr2) # Regular
                    ec += k*(-dr+2r)*(θr-1)/(θr+1)*0.5dr2/r # Due to ghost cell extrapolation in extrapolation
                    rhs[imx] -= Tf*k*(-0.5dr+2r) *dr2 /r/(θr+1) # Dirichlet BC in ghost cell extrap
                end
            elseif r_bulk && ir >= 2 && ir <= nr-1 
                ec +=  1.0dr2 + 0.5dr1*r1
                pc += -2.0dr2
                wc +=  1.0dr2 - 0.5dr1*r1
                rhs[imx] += 0
            else
                @warn "No r stencil!" iz ir
            end

            if n_front
                θz = ϕst[1] / (ϕst[1] - ϕst[3])
                if θz > dz
                    pc += -k*dz2*(θz+1)/θz
                    sc += k*dz2
                    rhs[imx] -= Tf*k*dz2/θz
                else
                    pc += -2k*dz2
                    sc += k*dz2*θz/(θz+1)
                    rhs[imx] -= 2Tf*k*dz2/(θz+1)
                end
            elseif s_front
                θz = ϕst[1] / (ϕst[1] - ϕst[5])
                if θz > dz
                    pc += -k*dz2*(θz+1)/θz
                    nc += k*dz2
                    rhs[imx] -= Tf*k*dz2/θz
                else
                    pc += -2k*dz2
                    nc += k*dz2*(2*θz)/(θz+1)
                    rhs[imx] -= 2Tf*k*dz2/(θz+1)
                end
            elseif z_bulk && iz >= 2 && iz <= nz-1
                sc +=  1.0dz2
                pc += -2.0dz2
                nc +=  1.0dz2
                rhs[imx] += 0
            else
                @warn "No z stencil!" iz ir
            end
            add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
            add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
            add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
            add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
            add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)
            rhs[imx] += - Q_ck 

            # push!(vals,                     1.0dz2); push!(cols, imx-nr); push!(rows, imx) # S cell
            # push!(vals,  1.0dr2 - 0.5dr1*r1       ); push!(cols, imx- 1); push!(rows, imx) # W cell
            # push!(vals, -2.0dr2            -2.0dz2); push!(cols, imx   ); push!(rows, imx) # P cell
            # push!(vals,  1.0dr2 + 0.5dr1*r1       ); push!(cols, imx+ 1); push!(rows, imx) # E cell
            # push!(vals,                     1.0dz2); push!(cols, imx+nr); push!(rows, imx) # N cell
        end
    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    sol = mat_lhs \ rhs
    T = reshape(sol, nr, nz)
end

function add_to_vcr!(vcr, dom, p_imx, shift, val)
    vals, cols, rows = vcr
    c_imx = p_imx + shift[1] + dom.nr*shift[2]
    push!(vals, val)
    push!(cols, c_imx)
    push!(rows, p_imx)
end


# solve_T = solve_T_original