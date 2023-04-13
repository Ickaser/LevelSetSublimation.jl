export make_artificial_params, make_default_params
export make_ϕ0
export params_nondim_setup

"""
    make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)

Return a ϕ0 with appropriate size for the passed Domain.

Parameter ϵ * sqrt(dom.rmax*dom.zmax) is added to ensure that interface is within domain.
Currently allowed setups:
- `:top`, `:flat` -- interface at zmax - ϵ
- `:rad`, `:cyl`  -- interface at rmax - ϵ
- `:box`          -- interface at both zmax-ϵ and rmax-ϵ
- `:cone`         -- interface a line decreasing in r
- `:ell_bub `     -- ellipse in center of vial, separated from boundaries
- `:circ `        -- circle at r=0, z=0 
- `:tinycirc `    -- circle at r=0, z=0, very small radius
"""
function make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)
    ϵ *= sqrt(dom.rmax * dom.zmax)
    if ϕtype == :top || ϕtype == :flat
        ϕ0 = [z - dom.zmax + ϵ for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :rad || ϕtype == :cyl
        ϕ0 = [r - dom.rmax + ϵ for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :box
        ϕ0 = [max(r-dom.rmax, z-dom.zmax) + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :ell_bub
        @unpack rmax, zmax = dom
        ϕ0 = [1.5* r^2 + 6*(z-0.5zmax)^2 - 1.0 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :circ
        ϕ0 = [sqrt(1.1r^2 / dom.rmax^2 + 1.1z^2 / (dom.zmax)^2 )   - 1.0 
                for r in dom.rgrid, z in dom.zgrid] .*sqrt(dom.rmax*dom.zmax)
    elseif ϕtype == :tinycirc
        ϕ0 = [dom.rmax * r^2 + dom.zmax * z^2 - 0.11 
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
    params_setup(cparams, controls)
    
Return a copied `params` dictionary and list of keys from `controls` 
to be updated at each time sampled point.
"""
function params_nondim_setup(cparams, controls)
    params = deepcopy(cparams)
    nondim_controls = deepcopy(controls)

    for pk in keys(cparams)
        try
            params[pk] = ustrip.(pbu[pk], cparams[pk])
        catch DomainError
            @error "Bad dimensions in passed parameter." pk params[pk] pbu[pk]   
        end
    end
    for mk in keys(controls)
        nondim_controls[mk] = ustrip.(pbu[mk], controls[mk])
    end


    arr_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])>1)]
    sc_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])==1)]
    for ki in sc_keys
        params[ki] = nondim_controls[ki]
    end
    for ki in arr_keys
        if length(controls[ki]) != length(controls[:t_samp])
            @error "Improper length of measurement vector. Must match time vector length." ki length(controls[ki]) length(controls[:t_samp])
        end
        params[ki] = nondim_controls[ki][1]
    end
    if length(arr_keys) == 0
        arr_keys = nothing
    end
    return params, arr_keys, nondim_controls
end

pbu = params_base_units = Dict{Symbol, Any}(
    :Kgl => u"W/m^2/K",
    :Kv => u"W/m^2/K",
    :Q_ic => u"W/m^3",
    :Q_gl_RF => u"W",
    :Q_ck => u"W/m^3",
    :k => u"W/m/K",

    :Tsh => u"K",
    :Tgl0 => u"K",
    :Tf0 => u"K",
    :Tgl => u"K",
    :Tf => u"K",

    :p_ch => u"Pa",

    :ρf => u"kg/m^3",
    :Cpf => u"J/kg/K",
    :ΔH => u"J/kg",
    :m_cp_gl => u"J/K",

    :ϵ => NoUnits,
    :l => u"m",
    :κ => u"m^2",
    :μ => u"Pa*s",
    :R => u"J/mol/K",
    :Mw => u"kg/mol",

    :t_samp => u"s",
)

"""
    make_default_params()

Return a dictionary of `cparams`, with values corresponding to SI units

In theory, gives physical values of parameters. Haven't actually done that, though.
"""
function make_default_params()
            
    # Very artificial parameters
    # Heat transfer
    Kgl = 1e2 *u"W/K/m^2" # 1/(Contact resistance 1e-4 + cylindrical resistance from outer wall 1e-3 to 1e-5)
    Kv = 20 * u"W/K/m^2"
    Q_ic = 0.3u"W/cm^3"
    Q_ck = 0.0u"W/m^3"
    k = k_sucrose
    m_cp_gl = 5u"g" * cp_gl # Half of a 10R vial's mass contributing; all of a 2R.

    # Mass transfer
    p_ch = 100u"mTorr" # 100 mTorr is about 13 Pa
    ϵ = 0.9 # 90% porosity
    l = 1e-6u"m" # ~size of a pore
    κ = 1e-10u"m^2" # ~size^2 of a pore
    R = 8.3145u"J/mol/K"
    Mw = .018u"kg/mol" #mol/kg
    μ = LevelSetSublimation.μ

    # Sublimation
    ΔH = LevelSetSublimation.ΔH
    # ΔHsub = 678.0 # u"cal/g"
    ρf = ρ_ice 
    Cpf = Cp_ice

    # Rw = 8.3145 / .018 # J/molK * mol/kg
    # calc_ρvap(T) = calc_psub(T)/Rw/T


    params = Dict{Symbol, Any}()
    @pack! params = Kgl, Kv, Q_ic, Q_ck, k, m_cp_gl, ΔH, ρf, p_ch, ϵ, l, κ, R, Mw, μ, Cpf
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
    Kgl = 1.0
    Kv = 1.0
    Q_ic = 1.0
    Q_ck = 0.0
    k = 1.0
    m_cp_gl = 5.0

    # Mass transfer
    p_ch = 5 # 100 mTorr is about 13 Pa
    ϵ = 0.9
    l = 1.0
    κ = 0.5
    R = 8.3145
    Mw = .018 #mol/kg
    μ = 1.0

    # Sublimation
    ΔH = 1.0
    # ΔHsub = 678.0 # u"cal/g"
    ρf = 100.0 
    Cpf = 10.0

    # Rw = 8.3145 / .018 # J/molK * mol/kg
    # calc_ρvap(T) = calc_psub(T)/Rw/T


    params = Dict{Symbol, Any}()
    @pack! params = Kgl, Kv, Q_ic, Q_ck, k, m_cp_gl, ΔH, ρf, p_ch, ϵ, l, κ, R, Mw, μ, Cpf
    return params
end

default_params = make_artificial_params()