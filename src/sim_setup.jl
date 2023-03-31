export make_artificial_params, make_decent_params
export make_ϕ0

"""
    make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)

Return a ϕ0 with appropriate size for the passed Domain.

Parameter ϵ is added to ensure that interface is within domain.
Currently allowed setups:
- `:top`, `:flat` -- interface at zmax - ϵ
- `:rad`, `:cyl`  -- interface at rmax - ϵ
- `:box`          -- interface at both zmax-ϵ and rmax-ϵ
- `:cone`         -- interface a line decreasing in r
- `:ell_bub `     -- ellipse in center of vial, separated from boundaries
- `:circ `        -- circle at r=0, z=0 
- `:circ_bub `    -- circle at r=0, z=0.5zmax, very small radius
"""
function make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)
    if ϕtype == :top || ϕtype == :flat
        ϕ0 = [z - dom.zmax + ϵ for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :rad || ϕtype == :cyl
        ϕ0 = [r - dom.rmax + ϵ for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :box
        ϕ0 = [max(r-dom.rmax, z-dom.zmax) + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :ell_bub
        @unpack rmax, zmax = dom
        ϕ0 = [1.5rmax * r^2 + 6zmax*(z-0.5zmax)^2 - 1.0 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :circ
        ϕ0 = [1.1dom.rmax * r^2 + 1.1dom.zmax * z^2 - 1.0 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :cone
        ϕ0 = [0.5r + z - 0.9dom.zmax + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    else
        @error "ArgumentError: Invalid ϕ0 kind to make_ϕ0" ϕtype
    end
    return ϕ0
end

"""
    make_decent_params()

Return a dictionary of `params`, with values corresponding to SI units

In theory, gives physical values of parameters. Haven't actually done that, though.
Also, currently broken.
"""
function make_decent_params()
    Q_gl = 2.0  # heat flux from glass
    Q_sh = 1.0  # heat flux from shelf
    Q_ic = 1.0  # volumetric heat in ice
    Q_ck = 0.0  # volumetric heat in cake
    k = 1.0     # cake thermal conductivity
    Tf = 250.0  # constant ice temperature
    ΔH = 10.0   # heat of sublimation
    # ΔHsub = 678.0 # u"cal/g"
    # ΔHsub = 2.837e9 # u"J/kg"
    ρf = 920    # density of ice
    params = Dict{Symbol, Any}()
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
    @pack! params = Q_gl, Q_sh, Q_ic, Q_ck, k, Tf, ΔH, ρf
    return params
end
        
"""
    make_artificial_params()

Return a dictionary of `T_params`, with artificially chosen values.

For convenience in testing code.
"""
function make_artificial_params()
            
    # Very artificial parameters
    # Heat transfer
    Q_gl = 1.0
    Q_sh = 1.0
    Q_ic = 1.0
    Q_ck = 0.0
    k = 1.0
    Tf = 250.0

    # Mass transfer
    p_sub = 5.0
    p_ch = 1.5
    ϵ = 0.9
    l = 1.0
    κ = 0.5
    R = 8.3145
    Mw = 18 
    μ = 1.0

    # Sublimation
    ΔH = 1.0
    # ΔHsub = 678.0 # u"cal/g"
    ρf = 100.0 

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


    params = Dict{Symbol, Any}()
    @pack! params = Q_gl, Q_sh, Q_ic, Q_ck, k, Tf, ΔH, ρf, p_sub, p_ch, ϵ, l, κ, R, Mw, μ
    return params
end

default_params = make_artificial_params()