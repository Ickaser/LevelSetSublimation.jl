# Construct sparse array with rows, cols, vals format, then construct sparse matrix
function solve_T(ϕ, dom, params)
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
        
function make_decent_params()
    Q_gl = 2.0  # heat flux from glass
    Q_sh = 1.0  # heat flux from shelf
    Q_ic = 1.0  # volumetric heat in ice
    Q_ck = 0.0  # volumetric heat in cake
    k = 1.0     # cake thermal conductivity
    Tf = 250.0  # constant ice temperature
    ΔH = 10.0   # heat of sublimation
    # ΔHsub = 678.0 # u"cal/g"
    ρf = 920    # density of ice
    T_params = Dict{Symbol, Any}()
    " Takes T as Kelvin, returns P in Pa"
    function calc_psub(T)
        ai = [-0.212144006e2,  0.273203819e2,  -0.610598130e1]
        bi = [0.333333333e-2,  0.120666667e1,  0.170333333e1]
        θ = T/275.16
        lnπ = sum(ai .* θ .^bi) / θ
        exp(lnπ)*611.657
    end
    Rw = 8.3145 / .018 # J/molK * mol/kg
    calc_ρvap(T) = calc_psub(T)/Rw/T # compute density of vapor
    @pack! T_params = Q_gl, Q_sh, Q_ic, Q_ck, k, Tf, ΔH, ρf
    return T_params
end
        
function make_artificial_params()
            
    # Very artificial parameters
    Q_gl = 1.0
    Q_sh = 1.0
    Q_ic = 1.0
    Q_ck = 0.0
    k = 1.0
    Tf = 250.0
    ΔH = 1.0
    # ΔHsub = 678.0 # u"cal/g"
    ρf = 100.0 
    T_params = Dict{Symbol, Any}()
    " Takes T as Kelvin, returns P in Pa"
    function calc_psub(T)
        ai = [-0.212144006e2,  0.273203819e2,  -0.610598130e1]
        bi = [0.333333333e-2,  0.120666667e1,  0.170333333e1]
        θ = T/275.16
        lnπ = sum(ai .* θ .^bi) / θ
        exp(lnπ)*611.657
    end
    Rw = 8.3145 / .018 # J/molK * mol/kg
    calc_ρvap(T) = calc_psub(T)/Rw/T
    @pack! T_params = Q_gl, Q_sh, Q_ic, Q_ck, k, Tf, ΔH, ρf
    return T_params
end

default_params = make_artificial_params()
