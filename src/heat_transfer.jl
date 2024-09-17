export solve_T

function calc_QpppRFd(params)
    base, tcp, tvps = params
    return 2π*tvps.f_RF*base.ε0*base.εpp_d*tcp.B_d*tvps.P_per_vial
end
function calc_QpppRFf(params)
    base, tcp, tvps = params
    return 2π*tvps.f_RF*base.ε0*base.εpp_f*tcp.B_f*tvps.P_per_vial
end
function calc_QpppRFvw(params)
    base, tcp, tvps = params
    return 2π*tvps.f_RF*base.ε0*base.εpp_vw*tcp.B_vw*tvps.P_per_vial
end

"""
    solve_T(u, Tf, dom::Domain, params)

Compute 2D axisymmetric T profile, returning ghost cell values, for given state `u`.

`u` contains ϕ and Tf; unpacked with `ϕ_T_from_u`.

Neumann boundary conditions on all rectangular boundaries; Dirichlet on zero-level set.

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
    @unpack Kvwf, kd = params[2]
    @unpack Kshf, Tsh = params[3]
    QRFd = calc_QpppRFd(params)
    # ϕ, Tf, Tvw = ϕ_T_from_u(u, dom)
    ϕ, Tvw = ϕ_T_from_u(u, dom)[[true, false, true]]

    # To prevent blowup, artificially add some  corner ice if none is present
    # This is by tampering with level set field, hopefully memory safe
    # if minimum(ϕ) > 0
    #     # @info "Solving heat equation without any ice, artificially introducing some"
    #     ϕ = copy(ϕ)
    #     ϕ[argmin(ϕ)] = - max(dr, dz)
    # end


    rows = Vector{Int}(undef, 0)
    cols = Vector{Int}(undef, 0)
    vals = similar(Tf, 0)
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
                    pc += -2kd*dr2/θr
                    rhs[imx] -= 2kd*Tf_loc*dr2/θr #+ BC1*(r1 - 2dr1)
                else # Front is within .05 cells of boundary
                    pc += -2kd*dr2/(θr+1)
                    rhs[imx] -= 2Tf_loc*kd*dr2/(θr+1) 
                end
            else
                # Using Neumann boundary to define ghost point: west T= east T - 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2kd*dr2
                ec +=  2kd*dr2
                rhs[imx] += 0 # r=0, so 1/r = NaN
            end
        elseif ir == nr
            # Robin BC: glass
            wϕ = ϕ[ir-1, iz]
            if wϕ < 0 # Front is within a cell of boundary
                θr = pϕ/(pϕ-wϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir-1]-Tf[ir])
                if θr >= θ_THRESH
                    pc += -2kd*dr2/θr - Kvwf*(r1 + 2dr1)
                    rhs[imx] -= 2Tf_loc*kd*dr2/θr + Kvwf*Tvw*(r1 + 2dr1)
                else 
                    # First, use Robin BC to define an east ghost cell
                    # THen, extrapolate across Stefan boundary using east ghost & Stefan
                    pc += (Kvwf*(-r1 -2dr1*θr) + kd*(r1*dr1 - 2dr2))/(θr+1)
                    rhs[imx] -= (Tf_loc*kd*(2dr2-r1*dr1) + Tvw*Kvwf*(r1+2dr1*θr))/(θr+1)
                end

            else
                # Using Robin boundary to define ghost point: east T= west T + 2BC1*dr
                # Also, the first derivative term is given exactly by BC1/r
                # p. 65, 66 of project notes
                pc += -2kd*dr2 - Kvwf*(r1 + 2dr1)
                wc +=  2kd*dr2
                rhs[imx] -= Kvwf*Tvw*(2dr1 + r1)
            end

        else # r direction bulk
            # Check for Stefan boundary
            eϕ = ϕ[ir+1, iz]
            wϕ = ϕ[ir-1, iz]
            if eϕ <= 0 && wϕ <= 0
                # Pretend in bulk, rather than treat two ghost cells. THis is a rare case
                ec +=  kd*(1.0dr2 + 0.5dr1*r1)
                pc += -kd*(2.0dr2)
                wc +=  kd*(1.0dr2 - 0.5dr1*r1)
                rhs[imx] += 0
            elseif eϕ <= 0 # East ghost cell, across front
                θr = pϕ / (pϕ - eϕ)
                Tf_loc = Tf[ir] + θr*(Tf[ir+1]-Tf[ir])
                if θr >= θ_THRESH
                    # Quadratic
                    pc += kd*(-2*dr2 -(1-θr)*dr1*r1)/θr
                    wc += kd*(-dr1*r1*θr + 2dr2)/(θr+1)# Regular + b gradient
                    rhs[imx] -= Tf_loc*kd*(2dr2 + dr1*r1)/θr/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -2kd*dr2 # Regular 
                    # pc += kd*(0.5dr1*r1+ dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    # wc += kd*(-0.5dr1*r1 + dr2) # Regular
                    # rhs[imx] -= Tf_loc*kd*(0.5*dr+r) *dr2 *r1/θr # Dirichlet BC in ghost cell extrap
                else
                    # Linear extrapolation, a cell out
                    pc += -2kd*dr2 # Regular
                    wc += kd*(-dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf_loc*kd*( dr1*r1 + 2dr2) /(θr+1) # Dirichlet BC in ghost cell extrap
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
                    pc += (kd*(-2*dr2+(1-θr)*dr1*r1))/θr
                    ec += (kd*( dr1*r1*θr + 2dr2))/(θr+1)# Regular + b gradient
                    rhs[imx] -= Tf_loc*(kd*(2dr2 - dr1*r1))/θr/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Linear extrapolation for ghost cell
                    # pc += -2kd*dr2 # Regular 
                    # pc += kd*(-0.5dr1*r1 + dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    # ec += kd*( 0.5dr1*r1 + dr2) # Regular
                    # rhs[imx] -= Tf_loc*kd*(-0.5dr1*r1+dr2)/θr # Dirichlet BC in ghost cell extrap
                    # if ir > 2 # Avoid singular at r=0
                    #     # Logarithmic extrapolation for ghost cell
                    #     pc += kd*(log(rΓ/rm)*0.5*r1*dr1 + log(rΓ*rm*r1*r1)*dr2)/log(r/rΓ)
                    #     ec += kd*(0.5r1*dr1 + dr2)
                    #     rhs[imx] -= Tf_loc*kd*(0.5dr1*r1*log(rm*r1) + log(r/rm)*dr2)/log(r/rΓ)
                    # else 
                    #     # Linear extrapolation for ghost cell
                    #     pc += -2kd*dr2 # Regular 
                    #     pc += kd*(-0.5dr1*r1 + dr2)*(θr-1)/θr # Due to ghost cell extrapolation
                    #     ec += kd*( 0.5dr1*r1 + dr2) # Regular
                    #     rhs[imx] -= Tf_loc*kd*(-0.5dr1*r1+dr2)/θr # Dirichlet BC in ghost cell extrap
                    # end
                else # Very small θ
                    # if ir > 2
                    #     # Logarithmic extrapolation
                    #     rp = rgrid[ir+1]
                    #     pc += -2kd*dr2 # Regular
                    #     ec += kd*(0.5dr1*r1*log(rm/rp) + log(rΓ*rΓ/rp/rm)*dr2)/log(rΓ/rp)
                    #     rhs[imx] -= Tf_loc*kd*(0.5dr1*r1*log(rp/rm) + log(rm/rp)*dr2)/log(rΓ/rp)
                    # else
                    #     # Linear extrapolation
                    #     pc += -2kd*dr2 # Regular
                    #     ec += kd*(dr1*r1 + 2θr*dr2)/(θr+1)  
                    #     rhs[imx] -= Tf_loc*kd*(-dr1*r1 + 2dr2)/(θr+1) # Dirichlet BC in ghost cell extrap
                    # end
                    # Linear extrapolation
                    pc += -2kd*dr2 # Regular
                    ec += kd*(dr1*r1 + 2θr*dr2)/(θr+1)  
                    rhs[imx] -= Tf_loc*kd*(-dr1*r1 + 2dr2)/(θr+1) # Dirichlet BC in ghost cell extrap
                    # Treat as constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
                end

            else # Bulk, not at front 
                ec +=  kd*(1.0dr2 + 0.5dr1*r1)
                pc += -kd*(2.0dr2)
                wc +=  kd*(1.0dr2 - 0.5dr1*r1)
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
                    pc += -2kd*dz2/θz - 2Kshf*dz1
                    rhs[imx] -= 2*Tf[ir]*kd*dz2/θz + 2Kshf*Tsh*dz1
                else
                    # First Robin ghost defined, then Stefan ghost extrap
                    pc += -2kd*dz2 - 2Kshf*dz1
                    nc +=  2kd*dz2
                    rhs[imx] -= 2Kshf*Tsh*dz1
                end
            else
                # Robin boundary
                pc += -kd*2dz2 - 2Kshf*dz1
                nc +=  kd*2dz2
                rhs[imx] -= 2Tsh*Kshf*dz1
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
                    pc += -2*kd*dz2/θz
                    rhs[imx] -= 2*Tf[ir]*kd*dz2/θz + BC4*2*dz1
                else
                    pc += -2kd*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*kd*dz2/(θz+1) 
                end
            else
                # Using Neumann boundary to define ghost point: east T= west T + 2BC1*dr
                # p. 65, 66 of project notes
                pc += -2kd*dr2
                sc +=  2kd*dr2
                rhs[imx] += -2kd*BC4*dz1
            end

        else # Bulk in z, still need to check for Stefan front
            nϕ = ϕ[ir, iz+1]
            sϕ = ϕ[ir, iz-1]
            if nϕ <= 0 && sϕ <= 0
                # Pretend in bulk, rather than treat two ghost cells. Rare case
                sc +=  kd*1.0dz2
                pc += -kd*2.0dz2
                nc +=  kd*1.0dz2
                rhs[imx] += 0
            elseif nϕ <= 0
                # stefan_debug = true
                θz = pϕ / (pϕ - nϕ)
                # println("θz=$θz, ir=$ir, iz = $iz, north")
                if θz >= θ_THRESH
                    # Quadratic
                    pc += (kd*(-2*dz2))/θz
                    sc += (kd*(2dz2))/(θz+1)# Regular + b gradient
                    rhs[imx] -= Tf[ir]*(kd*2dz2)/θz/(θz+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -kd*dz2*(θz+1)/θz
                    # sc += kd*dz2
                    # rhs[imx] -= Tf[ir]*kd*dz2/θz
                else
                    # Linear
                    pc += -2kd*dz2
                    sc += 2*kd*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*kd*dz2/(θz+1)
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
                    pc += (kd*(-2*dz2))/θz
                    nc += (kd*(2dz2))/(θz+1)# Regular + b gradient
                    rhs[imx] -= Tf[ir]*(kd*2dz2)/θz/(θz+1) # Dirichlet BC in ghost cell extrap
                    # Linear
                    # pc += -kd*dz2*(θz+1)/θz
                    # nc += kd*dz2
                    # rhs[imx] -= Tf[ir]*kd*dz2/θz
                else
                    # Linear
                    pc += -2kd*dz2
                    nc += 2*kd*θz*dz2/(θz+1)
                    rhs[imx] -= 2Tf[ir]*kd*dz2/(θz+1)
                    # Constant
                    # add_to_vcr!(vcr, dom, imx, ( 0, 0), 1)
                    # rhs[imx] = Tf[ir]
                    # continue
                end

            else # Bulk, no Stefan front
                sc +=  kd*1.0dz2
                pc += -kd*2.0dz2
                nc +=  kd*1.0dz2
                rhs[imx] += 0
            end
        end


        # Assign all computed stencil values into matrix
        pc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 0), pc)
        ec != 0 && add_to_vcr!(vcr, dom, imx, ( 1, 0), ec)
        wc != 0 && add_to_vcr!(vcr, dom, imx, (-1, 0), wc)
        nc != 0 && add_to_vcr!(vcr, dom, imx, ( 0, 1), nc)
        sc != 0 && add_to_vcr!(vcr, dom, imx, ( 0,-1), sc)
        rhs[imx] += - QRFd 


    end
    mat_lhs = sparse(rows, cols, vals, ntot, ntot)
    prob = LinearProblem(mat_lhs, rhs)
    sol = solve(prob, SparspakFactorization()).u 
    # sol = solve(prob, UMFPACKFactorization()).u 
    # sol = mat_lhs \ rhs
    T = reshape(sol, nr, nz)

    # if minimum(T) <= 0
    #     @info "negative temperatures" T
    # end
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
    Tfg = u[iTf(dom)]
    pseudosteady_Tf(u, dom, params, Tfg)
end
function pseudosteady_Tf(u, dom, params, Tf_g)
    ϕ, Tvw = ϕ_T_from_u(u, dom)[[true, false, true]]
    @unpack kf, ρf, Cpf = params[1]
    QRFf = calc_QpppRFf(params)
    dϕdx_all = dϕdx_all_WENO(ϕ, dom)
    has_ice = (compute_iceht_bottopcont(ϕ, dom)[1] .> 0)
    # Δξ = compute_iceht_bottopcont(ϕ, dom)[1]
    # @info "Tf solve" extrema(Δξ) extrema(Tf_g)
    if all(.~ has_ice) # If no ice present, skip nonlinear solve procedure
        return Tf_g
    end

    # function resid!(dTfdt, Tf)
    function resid!(dTfdt, Tf, unused_arg)
        if any(isnan.(Tf))
            @warn "NaN found" Tf
        end
        extrap_Tf_noice!(Tf, has_ice, dom)
        T = solve_T(u, Tf, dom, params)
        p = solve_p(u, Tf, T, dom, params)
        dTfdt_radial!(dTfdt, u, Tf, T, p, dϕdx_all, dom, params)

        nothing
    end

    if all(has_ice) # If all ice present, use all DOF

        prob = NonlinearProblem(resid!, Tf_g)
        sol = solve(prob, NewtonRaphson())
        # prob = SteadyStateProblem((du,u,unused,t)->resid!(du,u,unused), Tf_g)
        # sol = solve(prob, DynamicSS(Rosenbrock23()))

        # if sol.retcode == ReturnCode.MaxIters
        #     @info "maxit" sol.retcode sol.u
        #     prob_ss = SteadyStateProblem((du,u,unused,t)->resid!(du,u,unused), Tf_g)
        #     sol = solve(prob_ss, DynamicSS(Rosenbrock23()), maxiters=100)
        #     @info "SS iterations" sol.retcode sol.stats sol.u
        # end
        Tfs = sol.u
    else # If ice doesn't cover full radial extent, trim out those DOF
        Tf_trim = Tf_g[has_ice]
        
        function resid_lessdof!(dTfdt_trim, Tf_trim, unused_arg2)
        # resid! = function (dTfdt_trim, Tf_trim)
            dTfdt = zeros(eltype(dTfdt_trim), dom.nr)
            Tf = zeros(eltype(Tf_trim), dom.nr)
            Tf[has_ice] .= Tf_trim
            # Tf[.~has_ice] .= Tf_trim[1]
            extrap_Tf_noice!(Tf, has_ice, dom)
            resid!(dTfdt, Tf, unused_arg2)
            dTfdt_trim .= dTfdt[has_ice]
            nothing
        end

        prob = SteadyStateProblem((du,u,unused,t)->resid_lessdof!(du,u,unused), Tf_trim)
        sol = solve(prob, DynamicSS(Rosenbrock23()))
        # prob = NonlinearProblem(resid_lessdof!, Tf_trim)
        # sol = solve(prob, NewtonRaphson())

        Tfs = zeros(dom.nr)
        Tfs[has_ice] = sol.u
        extrap_Tf_noice!(Tfs, has_ice, dom)
    end

    return Tfs
end

extrap_quad(ir, ir1, ir2, ir3, Tf1, Tf2, Tf3) = Tf1*(ir-ir2)*(ir-ir3)/(ir1-ir2)/(ir1-ir3) + 
                Tf2*(ir-ir1)*(ir-ir3)/(ir2-ir1)/(ir2-ir3) +
                Tf3*(ir-ir1)*(ir-ir2)/(ir3-ir1)/(ir3-ir2)
extrap_lin(ir, ir1, ir2, Tf1, Tf2) = Tf1*(ir-ir2)/(ir1-ir2) + Tf2*(ir-ir1)/(ir2-ir1)

function right_extrap!(Tf, extrap_region, has_ice)
    left = extrap_region[begin]-1
    if checkbounds(Bool, has_ice, left-2) && has_ice[left-1] && has_ice[left-2]
        ref_pts = left .- (0:2)
        for ir in extrap_region
            Tf[ir] = extrap_quad(ir, ref_pts..., Tf[ref_pts]...)
        end
    elseif checkbounds(Bool, has_ice, left-1) && has_ice[left-1]
        ref_pts = left .- (0:1)
        for ir in extrap_region
            Tf[ir] = extrap_lin(ir, ref_pts..., Tf[ref_pts]...)
        end
    else
        Tf[extrap_region] .= Tf[left]
    end
end

function left_extrap!(Tf, extrap_region, has_ice)
    right = extrap_region[end]+1
    if checkbounds(Bool, has_ice, right+2) && has_ice[right+1] && has_ice[right+2]
        ref_pts = right .+ (0:2)
        for ir in extrap_region
            Tf[ir] = extrap_quad(ir, ref_pts..., Tf[ref_pts]...)
        end
    elseif checkbounds(Bool, has_ice, right+1) && has_ice[right+1]
        ref_pts = right .+ (0:1)
        for ir in extrap_region
            Tf[ir] = extrap_lin(ir, ref_pts..., Tf[ref_pts]...)
        end
    else
        Tf[extrap_region] .= Tf[right]
    end
end

function mid_extrap!(Tf, extrap_region, has_ice)
    if length(extrap_region) == 1
        loc = extrap_region[1]
        Tf[loc] = 0.5*(Tf[loc+1] + Tf[loc - 1])
    else
        n = length(extrap_region)
        leftsec = extrap_region[begin:n÷2]
        rightsec = extrap_region[n÷2+1:end]
        right_extrap!(Tf, leftsec, has_ice)
        left_extrap!(Tf, rightsec, has_ice)
    end
    if any(Tf[extrap_region] .< 0)
    end
end

function extrap_Tf_noice!(Tf, has_ice, dom)
    if all(.~ has_ice)
        return
    elseif sum(has_ice) == 1 # Only one ice cell: extrapolate to everywhere
        Tf[:] .= Tf[findfirst(has_ice)]
        return
    end

    edges = has_ice[1:end-1] .⊻ has_ice[2:end]

    if all(.~edges) # No gaps
        return # This case shouldn't actually be reached
    end

    # Treat first gap
    left = findfirst(edges)
    if ~has_ice[1] # First empty is at left edge, so no treatment necessary
        right = findfirst(has_ice)
        left_extrap!(Tf, 1:left, has_ice)
        # Extrapolation done
        edges[left] = false
        left = findfirst(edges)
    end
    while any(edges)
        edges[left] = false
        if all(.~edges) # Goes all the way to the boundary: right interpolate
            right_extrap!(Tf, left+1:dom.nr, has_ice)
            break
        end
        right = findfirst(edges)
        mid_extrap!(Tf, left+1:right, has_ice)
        edges[right] = false
        left = findfirst(edges) # If there are no edges left, this is fine
    end
    nothing
end

function compute_Qvwf(u, T, dom::Domain, params)
    @unpack Kvwf = params[2]
    ϕ, Tf, Tvw = ϕ_T_from_u(u, dom)
    # Heat flux from glass, at outer radius
    # zweights = compute_icegl_area_weights(ϕ, dom) # area for ice-glass
    zweights = fill(dom.dz, dom.nz)
    zweights[begin] = zweights[end] = dom.dz/2 # area for ice-glass + ice-cake
    Q_vwf = 2π*dom.rmax * Kvwf * sum(zweights .* ( Tvw .- T[end,:]))
end
