export make_artificial_params, make_decent_params
export make_ϕ0
export params_setup

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
- `:tinycirc `    -- circle at r=0, z=0, very small radius
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
function params_setup(cparams, controls)
    params = deepcopy(cparams)
    arr_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])>1)]
    sc_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])==1)]
    for ki in sc_keys
        params[ki] = controls[ki]
    end
    for ki in arr_keys
        if length(controls[ki]) != length(controls[:t_samp])
            @error "Improper length of measurement vector. Must match time vector length." ki length(controls[ki]) length(controls[:t_samp])
        end
        params[ki] = controls[ki][1]
    end
    if length(arr_keys) == 0
        arr_keys = nothing
    end
    return params, arr_keys
end


"""
    make_decent_params()

Return a dictionary of `params`, with values corresponding to SI units

In theory, gives physical values of parameters. Haven't actually done that, though.
Also, currently broken.
"""
function make_decent_params()
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