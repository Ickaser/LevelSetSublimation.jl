export solve_T, solve_T_original
# y


# Construct sparse array with rows, cols, vals format, then construct sparse matrix

function radial_Tf(u, dom, params) 
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Kgl, Kv, Q_ck, k, Tsh = params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
end 

"""
    solve_T(u, dom::Domain, params)

Compute 2D axisymmetric T profile, returning ghost cell values, for given state `u`.

`u` contains ϕ and Tf; unpacked with `ϕ_T_from_u`.

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
Neumann boundaries use a ghost point & BC to define ghost cell, then use same stencil as normal.
Coefficients computed in `gfm_extrap.ipynb`, using Sympy.  
(For higher order, see Gibou and Fedkiw, 2005, "A fourth order accurate discretization ... Laplace ... ")  
"""
function solve_T(u, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Kgl, Kv, Q_ck, k, Tsh = params
    ϕ, Tf, Tgl = ϕ_T_from_u(u, dom)
    # To prevent blowup, artificially add some  corner ice if none is present
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

    θ_thresh = dr / dom.rmax

    for iz in 1:nz, ir in 1:nr
        # Row position in matrix: r is small iteration, z is outer iteration
        imx = ir + (iz-1)*nr

        pϕ = ϕ[ir, iz]

        # Check if in frozen domain; if so, fix temperature
        if pϕ <= 0
            add_to_vcr!(vcr, dom, imx, (0, 0), 1)
            rhs[imx] = Tf
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
        # z = zgrid[iz]
        

        # Check if on boundaries or in bulk, set up BC or diffeq  
        # TODO 
        # Currently, boundaries are not checking for interface crossing.


        stefan_debug = false


        # R direction boundaries
        if ir == 1
            # Symmetry BC
            BC1 = 0
            eϕ = ϕ[ir+1, iz]
            # Check for Stefan front
            if eϕ < 0 # Front is within a cell of boundary
                # p. 65 of project notes
                θr = pϕ/(pϕ-eϕ)
                if θr >= θ_thresh
                    pc += -2k*dr2/θr
                    rhs[imx] -= 2k*Tf*dr2/θr #+ BC1*(r1 - 2dr1)
                else # Front is within dr/r cells of boundary
                    # Have an exact value given by BC + front
                    add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                    rhs[imx] = Tf - θr*BC1*dr
                    continue
                end
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
            # Robin BC: glass
            wϕ = ϕ[ir-1, iz]
            if wϕ < 0 # Front is within a cell of boundary
                # p. 119 of project notes
                θr = pϕ/(pϕ-wϕ)
                if θr >= θ_thresh
                    pc += -2k*dr2/θr - Kgl*(r1 + 2dr1*θr)
                    rhs[imx] -= 2Tf*k*dr2/θr + Kgl*Tgl*(r1 + 2dr1)
                else 
                    # First, use Robin BC to define an east ghost cell
                    # THen, extrapolate across Stefan boundary using east ghost & Stefan
                    pc += (Kgl*(-r1 -2dr1*θr) + k*(r1*dr1 - 2dr2))/(θr+1)
                    rhs[imx] -= (Tf*k*(2dr2-r1*dr1) + Tgl*Kgl*(r1+2dr1*θr))/(θr+1)
                end

            else
                # Using Robin boundary to define ghost point: east T= west T + 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2dr2 - Kgl*(r1 + 2dr1)
                wc +=  2dr2
                rhs[imx] -= Kgl*Tgl*(2dr1 + r1)
            end

        else # r direction bulk
            # Check for Stefan boundary
            eϕ = ϕ[ir+1, iz]
            wϕ = ϕ[ir-1, iz]
            if eϕ <= 0 && wϕ <= 0
                # Pretend in bulk, rather than treat two ghost cells. THis is a rare case
                ec +=  k*(1.0dr2 + 0.5dr1*r1)
                pc += -k*(2.0dr2)
                wc +=  k*(1.0dr2 - 0.5dr1*r1)
                rhs[imx] += 0
            elseif eϕ <= 0 # East ghost cell, across front
                θr = pϕ / (pϕ - eϕ)
                if θr >= θ_thresh
                    pc += -2k*dr2 # Regular 
                    pc += k*(0.5dr1*r1+ dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    wc += k*(-0.5dr1*r1 + dr2) # Regular
                    rhs[imx] -= Tf*k*(0.5*dr+r) *dr2 *r1/θr # Dirichlet BC in ghost cell extrap
                else
                    pc += -2k*dr2 # Regular
                    wc += k*(-dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf*k*( dr1*r1 + 2dr2) /(θr+1) # Dirichlet BC in ghost cell extrap
                end
            elseif wϕ <= 0 # West ghost cell across front
                θr = pϕ / (pϕ - wϕ)
                if θr >= θ_thresh # Regular magnitude θ
                    pc += -2k*dr2 # Regular 
                    pc += k*(-0.5dr1*r1 + dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    ec += k*( 0.5dr1*r1 + dr2) # Regular
                    rhs[imx] -= Tf*k*(-0.5dr1*r1+dr2)/θr # Dirichlet BC in ghost cell extrap
                else # Very small θ
                    pc += -2k*dr2 # Regular
                    ec += k*(dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf*k*(-dr1*r1 + 2dr2)/(θr+1) # Dirichlet BC in ghost cell extrap
                end

            else # Bulk, not at front 
                ec +=  k*(1.0dr2 + 0.5dr1*r1)
                pc += -k*(2.0dr2)
                wc +=  k*(1.0dr2 - 0.5dr1*r1)
                rhs[imx] += 0
            end
        end

        # z direction discretization

        # z direction boundaries
        if iz == 1
            # Robin BC: shelf
            nϕ = ϕ[ir, iz+1]
            # Check for Stefan front
            if nϕ < 0 # Front is within a cell of boundary
                # stefan_debug = true
                # p. 65 of project notes
                θz = pϕ/(pϕ-nϕ)
                if θz > θ_thresh
                    pc += -2k*dz2/θz - 2Kv*dz1
                    rhs[imx] -= 2*Tf*k*dz2/θz + 2Kv*Tsh*dz1
                else
                    # Have an exact value given by BC + front
                    add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                    rhs[imx] = (k*Tf + Kv*Tsh*θz*dz) / (k + Kv*θz*dz)
                    continue
                end
                # No way to treat other equations, so cut it here
            else
                # Using Neumann boundary to define ghost point: south T= north T - 2BC1*dr
                # p. 65, 66 of project notes
                # FINISH THIS: NOT CORRECT
                pc += -k*2dz2 - 2Kv*dz1
                nc +=  k*2dz2
                rhs[imx] -= 2Tsh*Kv*dz1
            end
        elseif iz == nz
            # Adiabatic BC
            BC4 = 0
            sϕ = ϕ[ir, iz-1]
            if sϕ < 0 # Front is within a cell of boundary
                # stefan_debug = true
                # p. 65 of project notes
                θz = pϕ/(pϕ-sϕ)
                if θz > θ_thresh
                    pc += -2*k*dz2/θz
                    rhs[imx] -= 2*Tf*k*dz2/θz + BC4*2*dz1
                else
                    # # Have an exact value given by BC + front
                    # # No way to treat other equations, so add to matrix and stop here
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1) # P cell
                    # rhs[imx] = Tf + θz*BC4*dz
                    # continue
                    pc += -2k*dz2/(θz+1)
                    rhs[imx] -= 2Tf*k*dz2/(θz+1) 
                end
            else
                # Using Neumann boundary to define ghost point: east T= west T + 2BC1*dr
                # p. 65, 66 of project notes
                pc += -2k*dr2
                sc +=  2k*dr2
                rhs[imx] += -2k*BC4*dz1
            end

        else # Bulk in z, still need to check for Stefan front
            nϕ = ϕ[ir, iz+1]
            sϕ = ϕ[ir, iz-1]
            if nϕ <= 0 && sϕ <= 0
                # Pretend in bulk, rather than treat two ghost cells. Rare case
                sc +=  k*1.0dz2
                pc += -k*2.0dz2
                nc +=  k*1.0dz2
                rhs[imx] += 0
            elseif nϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - nϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, north")
                if θz >= θ_thresh
                    pc += -k*dz2*(θz+1)/θz
                    sc += k*dz2
                    rhs[imx] -= Tf*k*dz2/θz
                else
                    pc += -2k*dz2
                    sc += 2*k*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf*k*dz2/(θz+1)
                end
            elseif sϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - sϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, south")
                if θz >= θ_thresh
                    pc += -k*dz2*(θz+1)/θz
                    nc += k*dz2
                    rhs[imx] -= Tf*k*dz2/θz
                else
                    pc += -2k*dz2
                    nc += 2*k*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf*k*dz2/(θz+1)
                end

            else # Bulk, no Stefan front
                sc +=  k*1.0dz2
                pc += -k*2.0dz2
                nc +=  k*1.0dz2
                rhs[imx] += 0
            end
        end

        # Stefan boundary debug
        # if stefan_debug
        #     println("ir=$ir, iz=$iz, pc=$pc, nc=$nc, sc=$sc,  ec=$ec, wc=$wc, rhs=$(rhs[imx])") 
        # end

        # Debug boundaries
        # if iz==1
        #     # println("ir=$ir, iz=$iz, pc=$pc, nc=$nc, sc=$sc,  ec=$ec, wc=$wc, rhs=$(rhs[imx])") 
        # elseif iz == nz
        #     println("ir=$ir, iz=$iz, pc=$pc, nc=$nc, sc=$sc,  ec=$ec, wc=$wc, rhs=$(rhs[imx])") 
        # end

        # Assign all computed stencil values into matrix
        pc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
        ec != 0 && add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
        wc != 0 && add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
        nc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
        sc != 0 && add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)
        rhs[imx] += - Q_ck 

        # push!(vals,                     1.0dz2); push!(cols, imx-nr); push!(rows, imx) # S cell
        # push!(vals,  1.0dr2 - 0.5dr1*r1       ); push!(cols, imx- 1); push!(rows, imx) # W cell
        # push!(vals, -2.0dr2            -2.0dz2); push!(cols, imx   ); push!(rows, imx) # P cell
        # push!(vals,  1.0dr2 + 0.5dr1*r1       ); push!(cols, imx+ 1); push!(rows, imx) # E cell
        # push!(vals,                     1.0dz2); push!(cols, imx+nr); push!(rows, imx) # N cell
    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    sol = mat_lhs \ rhs
    T = reshape(sol, nr, nz)
end

function add_to_vcr!(vcr, dom, p_imx, shift, val)
    vals, cols, rows = vcr
    c_imx = p_imx + shift[1] + dom.nr*shift[2]
    # println("p_imx = $p_imx, c_imx = $c_imx, shift=$shift")
    push!(vals, val)
    push!(cols, c_imx)
    push!(rows, p_imx)
end


# solve_T = solve_T_original