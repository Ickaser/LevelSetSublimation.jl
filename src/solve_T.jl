export solve_T


# Construct sparse array with rows, cols, vals format, then construct sparse matrix
"""
    solve_T(ϕ, dom::Domain, params)
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
        