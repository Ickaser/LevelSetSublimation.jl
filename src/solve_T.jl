export solve_T

# function radial_Tf(u, dom, params) 
#     @unpack dr, dz, dr1, dz1, dr2, dz2, 
#             rgrid, zgrid, nr, nz, ntot = dom
#     @unpack Kw, Kv, Q_ck, k, Tsh = params
#     ϕ, Tf, Tw = ϕ_T_from_u(u, dom)
# end 

"""
    solve_T(u, Tf, dom::Domain, params)

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
function solve_T(u, Tf, dom::Domain, params)
    @unpack dr, dz, dr1, dz1, dr2, dz2, 
            rgrid, zgrid, nr, nz, ntot = dom
    @unpack Kw, Kv, Q_ck, k, Tsh = params
    # ϕ, Tf, Tw = ϕ_T_from_u(u, dom)
    ϕ, Tw = ϕ_T_from_u(u, dom)[[true, false, true]]

    # To prevent blowup, artificially add some  corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    # if minimum(ϕ) > 0
    #     # @info "Solving heat equation without any ice, artificially introducing some"
    #     ϕ = copy(ϕ)
    #     ϕ[argmin(ϕ)] = - max(dr, dz)
    # end


    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = Vector{eltype(Tf)}(undef, 0)

    vcr = (vals, cols, rows)
    rhs = similar(Tf, ntot)
    rhs .= 0


    for iz in 1:nz, ir in 1:nr
        # Row position in matrix: r is small iteration, z is outer iteration
        imx = ir + (iz-1)*nr

        pϕ = ϕ[ir, iz]
        # Check if in frozen domain; if so, fix temperature
        if pϕ <= 0
            add_to_vcr!(vcr, dom, imx, (0, 0), 1)
            rhs[imx] = Tf[ir]
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
                θr = pϕ/(pϕ-eϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                if θr >= θ_THRESH
                    pc += -2k*dr2/θr
                    rhs[imx] -= 2k*Tf_loc*dr2/θr #+ BC1*(r1 - 2dr1)
                else # Front is within .05 cells of boundary
                    pc += -2k*dr2/(θr+1)
                    rhs[imx] -= 2Tf_loc*k*dr2/(θr+1) 
                end
            else
                # Using Neumann boundary to define ghost point: west T= east T - 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2k*dr2
                ec +=  2k*dr2
                rhs[imx] += 0 # r=0, so 1/r = NaN
            end
        elseif ir == nr
            # Robin BC: glass
            wϕ = ϕ[ir-1, iz]
            if wϕ < 0 # Front is within a cell of boundary
                θr = pϕ/(pϕ-wϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                if θr >= θ_THRESH
                    pc += -2k*dr2/θr - Kw*(r1 + 2dr1)
                    rhs[imx] -= 2Tf_loc*k*dr2/θr + Kw*Tw*(r1 + 2dr1)
                else 
                    # First, use Robin BC to define an east ghost cell
                    # THen, extrapolate across Stefan boundary using east ghost & Stefan
                    pc += (Kw*(-r1 -2dr1*θr) + k*(r1*dr1 - 2dr2))/(θr+1)
                    rhs[imx] -= (Tf_loc*k*(2dr2-r1*dr1) + Tw*Kw*(r1+2dr1*θr))/(θr+1)
                end

            else
                # Using Robin boundary to define ghost point: east T= west T + 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2k*dr2 - Kw*(r1 + 2dr1)
                wc +=  2k*dr2
                rhs[imx] -= Kw*Tw*(2dr1 + r1)
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
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                if θr >= θ_THRESH
                    # Quadratic
                    pc += k*(-2*dr2 -(1-θr)*dr1*r1)/θr
                    wc += k*(-dr1*r1*θr + 2dr2)/(θr+1)# Regular + b gradient
                    rhs[imx] -= Tf_loc*k*(2dr2 + dr1*r1)/θr/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -2k*dr2 # Regular 
                    # pc += k*(0.5dr1*r1+ dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    # wc += k*(-0.5dr1*r1 + dr2) # Regular
                    # rhs[imx] -= Tf_loc*k*(0.5*dr+r) *dr2 *r1/θr # Dirichlet BC in ghost cell extrap
                else
                    # Linear extrapolation, a cell out
                    pc += -2k*dr2 # Regular
                    wc += k*(-dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf_loc*k*( dr1*r1 + 2dr2) /(θr+1) # Dirichlet BC in ghost cell extrap
                    # Constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
                end
            elseif wϕ <= 0 # West ghost cell across front
                θr = pϕ / (pϕ - wϕ)
                # rΓ = r - θr*dr
                # rm = rgrid[ir-1]
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                if θr >= θ_THRESH # Regular magnitude θ
                    # Quadratic ghost cell extrapolation
                    pc += (k*(-2*dr2+(1-θr)*dr1*r1))/θr
                    ec += (k*( dr1*r1*θr + 2dr2))/(θr+1)# Regular + b gradient
                    rhs[imx] -= Tf_loc*(k*(2dr2 - dr1*r1))/θr/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Linear extrapolation for ghost cell
                    # pc += -2k*dr2 # Regular 
                    # pc += k*(-0.5dr1*r1 + dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    # ec += k*( 0.5dr1*r1 + dr2) # Regular
                    # rhs[imx] -= Tf_loc*k*(-0.5dr1*r1+dr2)/θr # Dirichlet BC in ghost cell extrap
                    # if ir > 2 # Avoid singular at r=0
                    #     # Logarithmic extrapolation for ghost cell
                    #     pc += k*(log(rΓ/rm)*0.5*r1*dr1 + log(rΓ*rm*r1*r1)*dr2)/log(r/rΓ)
                    #     ec += k*(0.5r1*dr1 + dr2)
                    #     rhs[imx] -= Tf_loc*k*(0.5dr1*r1*log(rm*r1) + log(r/rm)*dr2)/log(r/rΓ)
                    # else 
                    #     # Linear extrapolation for ghost cell
                    #     pc += -2k*dr2 # Regular 
                    #     pc += k*(-0.5dr1*r1 + dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    #     ec += k*( 0.5dr1*r1 + dr2) # Regular
                    #     rhs[imx] -= Tf_loc*k*(-0.5dr1*r1+dr2)/θr # Dirichlet BC in ghost cell extrap
                    # end
                else # Very small θ
                    # if ir > 2
                    #     # Logarithmic extrapolation
                    #     rp = rgrid[ir+1]
                    #     pc += -2k*dr2 # Regular
                    #     ec += k*(0.5dr1*r1*log(rm/rp) + log(rΓ*rΓ/rp/rm)*dr2)/log(rΓ/rp)
                    #     rhs[imx] -= Tf_loc*k*(0.5dr1*r1*log(rp/rm) + log(rm/rp)*dr2)/log(rΓ/rp)
                    # else
                    #     # Linear extrapolation
                    #     pc += -2k*dr2 # Regular
                    #     ec += k*(dr1*r1 + 2θr*dr2)/(θr+1)  
                    #     rhs[imx] -= Tf_loc*k*(-dr1*r1 + 2dr2)/(θr+1) # Dirichlet BC in ghost cell extrap
                    # end
                    # Linear extrapolation
                    pc += -2k*dr2 # Regular
                    ec += k*(dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf_loc*k*(-dr1*r1 + 2dr2)/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Treat as constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
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
                if θz > θ_THRESH
                    pc += -2k*dz2/θz - 2Kv*dz1
                    rhs[imx] -= 2*Tf[ir]*k*dz2/θz + 2Kv*Tsh*dz1
                else
                    # First Robin ghost defined, then Stefan ghost extrap
                    pc += -2k*dz2 - 2Kv*dz1
                    nc +=  2k*dz2
                    rhs[imx] -= 2Kv*Tsh*dz1
                end
            else
                # Robin boundary
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
                if θz > θ_THRESH
                    pc += -2*k*dz2/θz
                    rhs[imx] -= 2*Tf[ir]*k*dz2/θz + BC4*2*dz1
                else
                    pc += -2k*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*k*dz2/(θz+1) 
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
                if θz >= θ_THRESH
                    # Quadratic
                    pc += (k*(-2*dz2))/θz
                    sc += (k*(2dz2))/(θz+1)# Regular + b gradient
                    rhs[imx] -= Tf[ir]*(k*2dz2)/θz/(θz+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -k*dz2*(θz+1)/θz
                    # sc += k*dz2
                    # rhs[imx] -= Tf[ir]*k*dz2/θz
                else
                    # Linear
                    pc += -2k*dz2
                    sc += 2*k*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*k*dz2/(θz+1)
                    # Constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
                end
            elseif sϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - sϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, south")
                if θz >= θ_THRESH
                    # Quadratic
                    pc += (k*(-2*dz2))/θz
                    nc += (k*(2dz2))/(θz+1)# Regular + b gradient
                    rhs[imx] -= Tf[ir]*(k*2dz2)/θz/(θz+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -k*dz2*(θz+1)/θz
                    # nc += k*dz2
                    # rhs[imx] -= Tf[ir]*k*dz2/θz
                else
                    # Linear
                    pc += -2k*dz2
                    nc += 2*k*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*k*dz2/(θz+1)
                    # Constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
                end

            else # Bulk, no Stefan front
                sc +=  k*1.0dz2
                pc += -k*2.0dz2
                nc +=  k*1.0dz2
                rhs[imx] += 0
            end
        end


        # Assign all computed stencil values into matrix
        pc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
        ec != 0 && add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
        wc != 0 && add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
        nc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
        sc != 0 && add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)
        rhs[imx] += - Q_ck 


    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    prob = LinearProblem(mat_lhs, rhs)
    sol = solve(prob, SparspakFactorization()).u 
    # sol = solve(prob, UMFPACKFactorization()).u 
    # sol = mat_lhs \ rhs
    T = reshape(sol, nr, nz)

    if minimum(T) <= 0
        @info "negative temperatures" T
    end
    return T
end

function add_to_vcr!(vcr, dom, p_imx, shift, val)
    vals, cols, rows = vcr
    c_imx = p_imx + shift[1] + dom.nr*shift[2]
    push!(vals, val)
    push!(cols, c_imx)
    push!(rows, p_imx)
end

function pseudosteady_Tf(u, dom, params)
    # ϕv, Tfv, Twv = ϕ_T_from_u(u, dom)#[[true, false, true]]
    Tfg = ϕ_T_from_u(u, dom)[2]
    pseudosteady_Tf(u, dom, params, Tfg)
end
function pseudosteady_Tf(u, dom, params, Tf_g)
    ϕ, Tw = ϕ_T_from_u(u, dom)[[true, false, true]]
    @unpack ρf, Cpf, kf, Q_ic = params
    @unpack kf, ρf, Cpf = params
    dϕdx_all = dϕdx_all_WENO(ϕ, dom)
    has_ice = (compute_iceht_bottopcont(ϕ, dom)[1] .> 0)
    if all(.~ has_ice) # If no ice present, skip nonlinear solve procedure
        return Tf_g
    end
    # if any(isnan.(Tf_g))
    #     @warn "NaN guess" Tf_g
    #     Tf_g[isnan.(Tf_g)] .= 245.0
    # end


    if all(has_ice) # IF all ice present, use all DOF
        resid! = function (dTfdt, Tf)
            if any(isnan.(Tf))
                @warn "NaN found" Tf
            end
            extrap_Tf_noice!(Tf, has_ice, dom)
            if any(clamp.(Tf[has_ice], 200, 300) .!= Tf[has_ice])
                if sum(has_ice) > 0.75*dom.nr
                    els = findall(.~has_ice)
                else
                    els = has_ice
                end
                if typeof(Tf[1]) <: AbstractFloat
                    @info "Crazy Tf" Tf[has_ice] els
                else
                    @info "Crazy Tf" [Tfi.value for Tfi in Tf][has_ice] els
                end
                clamp!(Tf[has_ice], 200, 300)
            end
            T = solve_T(u, Tf, dom, params)
            p = solve_p(u, Tf, T, dom, params)
            dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom, params)
        end
        sol = nlsolve(resid!, Tf_g, autodiff=:forward, ftol=1e-10)
        Tfs = sol.zero
    else # If ice doesn't cover full radial extent, trim out those DOF
        Tf_trim = Tf_g[has_ice]
        resid! = function (dTfdt_trim, Tf_trim)
            dTfdt = zeros(eltype(dTfdt_trim), dom.nr)
            Tf = zeros(eltype(Tf_trim), dom.nr)
            Tf[has_ice] .= Tf_trim
        # function resid!(dTfdt, Tf, ssparams, t)
            # Tfv .= Tf
            if any(isnan.(Tf))
                @warn "NaN found" Tf
            end
            extrap_Tf_noice!(Tf, has_ice, dom)
            if any(clamp.(Tf[has_ice], 200, 400) .!= Tf[has_ice])
                if sum(has_ice) > 0.75*dom.nr
                    els = findall(.~has_ice)
                else
                    els = has_ice
                end
                if typeof(Tf[1]) <: AbstractFloat
                    @info "Crazy Tf" Tf[has_ice] els
                else
                    @info "Crazy Tf" [Tfi.value for Tfi in Tf][has_ice] els
                end
                clamp!(Tf[has_ice], 200, 300)
            end
            # @info "resid"
            # @info "resid: $(norm(dTfdt, 1)), $(norm(dTfdt, Inf))"
            # T_cache .= solve_T(u, Tf, dom, params)
            # p_cache .= solve_p(u, Tf, T_cache, dom, params, p_cache)
            # dTfdt_radial!(dTfdt, u, Tf, T_cache, p_cache, dϕdx_all, dom, params)
            # Tf_extrap = copy(Tf)
            T = solve_T(u, Tf, dom, params)
            # Tm = fill(250.15, size(dom))
            p = solve_p(u, Tf, T, dom, params)
            # p = solve_p(u, Tf, Tm, dom, params)
            dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom, params)
            # dTfdt[no_ice] = Tfv[no_ice] .- Tf[no_ice]
            # @info "resid" extrema(dTfdt) extrema(Tf) extrema(Tf_extrap)
            dTfdt_trim .= dTfdt[has_ice]
            nothing
        end
        sol = nlsolve(resid!, Tf_trim, autodiff=:forward, ftol=1e-10)
        Tfs = zeros(dom.nr)
        Tfs[has_ice] = sol.zero
        extrap_Tf_noice!(Tfs, has_ice, dom)
    end

    return Tfs
end

function pseudosteady_Tf_T_p(u, dom, params; abstol=1e-2)
    Tf0 = ϕ_T_from_u(u, dom)[2]
    T0 = solve_T(u, Tf0, dom, params)
    p0 = solve_p(u, Tf0, T0, dom, params)
    return pseudosteady_Tf_T_p(u, dom, params, Tf0, p0; abstol=abstol)
end

function pseudosteady_Tf_T_p(u, dom, params, Tfg, pg; abstol=1e-2)
    ϕ, Tf, Tw = ϕ_T_from_u_view(u, dom)
    @unpack kf, ρf, Cpf = params

    dϕdx_all = dϕdx_all_WENO(ϕ, dom)

    Tf .= Tfg
    T0 = solve_T(u, Tf, dom, params)
    p0 = solve_p(u, Tf, T0, dom, params, pg)
    Δξ, bot_contact, top_contact = compute_iceht_bottopcont(ϕ, dom)
    α = kf/ρf/Cpf
    CFL = 0.4
    dt = CFL / (α/dom.dr^2) #* max(0.01, minimum(Δξ)/dom.zmax)
    nt = max(3*dom.nr, ceil(Int, 100.0/dt))
    has_ice = (Δξ .> 0)

    dTfdt = zeros(dom.nr)

    for it in 1:nt
        # @info "step: it=$it, ext=$(extrema(dTfdt))"
        dTfdt_radial!(dTfdt, u, Tf, T0, p0, dϕdx_all, dom, params)
        @. Tf += dTfdt*dt
        extrap_Tf_noice!(Tf, has_ice, dom)
        T0 .= solve_T(u, Tf, dom, params)
        p0 .= solve_p(u, Tf, T0, dom, params, p0)
        if norm(dTfdt, Inf) < abstol
            @info "num steps to pseudosteady" it
            # break
            return Tf, T0, p0
        end
    end




    @info "max steps reached" nt dTfdt
    return Tf, T0, p0

end

function extrap_Tf_noice!(Tf, has_ice, dom)
    if all(.~ has_ice)
        return
    elseif sum(has_ice) == 1 # Only one ice cell: extrapolate to everywhere
        Tf[:] .= Tf[findfirst(has_ice)]
        return
    end

    first = findfirst(has_ice)
    if first > 1 && first+1<dom.nr && has_ice[first+1] && has_ice[first+2]
        # Enough points to do a quadratic extrapolation, so do that.
        # Is fine because we only really need one grid point of extrapolation
        ir1 = findfirst(has_ice)
        ir2 = ir1 + 1
        ir3 = ir1 + 2
        Tf1 = Tf[ir1]
        Tf2 = Tf[ir2]
        Tf3 = Tf[ir3]
        left_Textrap_quad(ir) = Tf1*(ir-ir2)*(ir-ir3)/(ir1-ir2)/(ir1-ir3) + 
                           Tf2*(ir-ir1)*(ir-ir3)/(ir2-ir1)/(ir2-ir3) +
                           Tf3*(ir-ir1)*(ir-ir2)/(ir3-ir1)/(ir3-ir2)
        for ir in 1:ir1-1
            Tf[ir] = left_Textrap_quad(ir)
        end
    elseif first > 1 && first < dom.nr && has_ice[first+1]
        ir1 = findfirst(has_ice)
        ir2 = ir1 + 1
        Tf1 = Tf[ir1]
        Tf2 = Tf[ir2]
        # Assuming uniform grid, can work in indices rather than space
        # Build a linear extrapolation
        left_Textrap(ir) = (Tf2-Tf1)/(ir2-ir1) * (ir-ir1) + Tf1
        for ir in 1:ir1-1
            Tf[ir] = left_Textrap(ir)
        end
    elseif first > 1 # No neighboring ice, so constant extrapolation
        Tf[1:findfirst(has_ice)-1] .= Tf[findfirst(has_ice)]
    end
    last = findlast(has_ice)
    if last < dom.nr && last-1>1 && has_ice[last-1] && has_ice[last-2]
        # Enough points to do a quadratic extrapolation, so do that.
        # Is fine because we only really need one grid point of extrapolation
        ir1 = findlast(has_ice)
        ir2 = ir1 - 1
        ir3 = ir1 - 2
        Tf1 = Tf[ir1]
        Tf2 = Tf[ir2]
        Tf3 = Tf[ir3]
        right_Textrap_quad(ir) = Tf1*(ir-ir2)*(ir-ir3)/(ir1-ir2)/(ir1-ir3) + 
                            Tf2*(ir-ir1)*(ir-ir3)/(ir2-ir1)/(ir2-ir3) +
                            Tf3*(ir-ir1)*(ir-ir2)/(ir3-ir1)/(ir3-ir2)
        for ir in ir1+1:dom.nr
            Tf[ir] = right_Textrap_quad(ir)
        end
        # typeof(Tf1) <: AbstractFloat && @info "R extrap" Tf[ir3:ir1+2]
    elseif last < dom.nr && last > 1 && has_ice[last-1]
        # Build a linear extrapolation
        ir1 = findlast(has_ice)
        ir2 = ir1 - 1
        Tf1 = Tf[ir1]
        Tf2 = Tf[ir2]
        # Assuming uniform grid, can work in indices rather than space
        right_Textrap(ir) = (Tf2-Tf1)/(ir2-ir1) * (ir-ir1) + Tf1
        for ir in ir1+1:dom.nr
            Tf[ir] = right_Textrap(ir)
        end
    elseif last < dom.nr # No neighboring ice, so constant extrapolation
        Tf[findlast(has_ice)+1:end] .= Tf[findlast(has_ice)]
    end

    edges = has_ice[1:end-1] .⊻ has_ice[2:end]
    if sum(edges) > 2
        gaps = findall(edges)[2:end-1] 
        for g in gaps
            if g-1 ∉ gaps && g+1 ∉ gaps
                Tf[g] = (Tf[g-1] + Tf[g+1]) / 2
            else
                @warn "unhandled: large gaps in ice"
            end
        end
    end
end
