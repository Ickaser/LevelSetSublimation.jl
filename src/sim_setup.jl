export make_default_params
export make_ϕ0, make_ramp
export params_nondim_setup, make_u0_ndim
export ϕ_T_from_u, ϕ_T_from_u_view

"""
    make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)

Return a ϕ0 with appropriate size for the passed Domain.

Parameter ϵ * sqrt(dom.rmax*dom.zmax) is added to ensure that interface is within domain.
Currently allowed setups:
- `:top`, `:flat` -- interface at zmax*(1 - ϵ)
- `:rad`, `:cyl`  -- interface at rmax*(1 - ϵ)
- `:box`          -- interface at both zmax*(1-ϵ) and rmax*(1-ϵ)
- `:cone`         -- interface a line decreasing in r
- `:midflat`      -- interface at zmax*0.5
- `:ell_bub `     -- ellipse in center of vial, separated from boundaries
- `:circ `        -- circle at r=0, z=0 
- `:tinycirc `    -- circle at r=0, z=0, very small radius
"""
function make_ϕ0(ϕtype::Symbol, dom::Domain; ϵ=1e-4)
    # ϵ *= sqrt(dom.rmax * dom.zmax)
    if ϕtype == :top || ϕtype == :flat
        ϕ0 = [z - dom.zmax*(1-ϵ) for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :midflat
        ϕ0 = [z - dom.zmax*0.5 for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :rad || ϕtype == :cyl
        ϕ0 = [r - dom.rmax*(1-ϵ) for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :box
        ϕ0 = [max(r-dom.rmax*(1-ϵ), z-dom.zmax*(1-ϵ))  
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :circ
        ϕ0 = [1.1*sqrt(r^2 / dom.rmax^2 + z^2 / (dom.zmax)^2 ) - 1.0
                for r in dom.rgrid, z in dom.zgrid] .*sqrt(dom.rmax*dom.zmax)
    elseif ϕtype == :circ_end
        @warn "The more arcane starting shapes are not necessarily well-tested."
        ϕ0 = [1.1*sqrt(r^2 / dom.rmax^2 + (z+0.75dom.zmax)^2 / (dom.zmax)^2 ) - 1.0
                for r in dom.rgrid, z in dom.zgrid] .*sqrt(dom.rmax*dom.zmax) 
    elseif ϕtype == :ell_bub
        @warn "The more arcane starting shapes are not well-tested."
        @unpack rmax, zmax = dom
        ϕ0 = [1.5* r^2 + 6*(z-0.5zmax)^2 - 1.0 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :tinycirc
        @warn "The more arcane starting shapes are not well-tested."
        ϕ0 = [sqrt(10r^2 / dom.rmax^2 + 10z^2 / (dom.zmax)^2 )   - 1.0 
                for r in dom.rgrid, z in dom.zgrid] .*sqrt(dom.rmax*dom.zmax)
    elseif ϕtype == :cone
        @warn "The more arcane starting shapes are not well-tested."
        ϕ0 = [0.5r + z - 0.9dom.zmax + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    else
        @error "ArgumentError: Invalid ϕ0 kind to make_ϕ0" ϕtype
    end
    return ϕ0
end

"""
    params_nondim_setup(cparams, controls)
    
Return a copied `params` dictionary and list of keys from `controls` 
to be updated at each time sampled point.
"""
function params_nondim_setup(cparams, controls)
    params = deepcopy(cparams)
    nondim_controls = deepcopy(controls)

    for pk in keys(cparams)
        try
            params[pk] = ustrip.(PBD[pk], cparams[pk])
        catch DomainError
            @error "Bad dimensions in passed parameter." pk params[pk] PBD[pk]   
        end
    end
    for mk in keys(controls)
        # nondim_controls[mk] = ustrip.(PBD[mk], controls[mk])
        nondim_controls[mk] = nondim_controlvar(mk, controls[mk])
        params[mk] = nondim_controls[mk](0)
    end
    # if length(nondim_controls[:t_samp]) == 1
    #     nondim_controls[:t_samp] = [0.0]
    # end

    # arr_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])>1)]
    # sc_keys = [ki for ki in keys(controls) if (ki!=:t_samp && length(controls[ki])==1)]
    # for ki in sc_keys
    #     params[ki] = nondim_controls[ki]
    # end
    # for ki in arr_keys
    #     # if length(controls[ki]) != length(controls[:t_samp])
    #     #     @error "Improper length of measurement vector. Must match time vector length." ki length(controls[ki]) length(controls[:t_samp])
    #     # end
    #     params[ki] = nondim_controls[ki](0)
    # end
    # if length(arr_keys) == 0
    #     arr_keys = nothing
    # end
    # return params, arr_keys, nondim_controls

    return params, nondim_controls
end

function nondim_controlvar(varname::Symbol, control_dim::RampedVariable)
    if dimension(control_dim(0u"s")) != dimension(PBD[varname])
        @error "Bad units on ramped variable." varname control_dim PBD[varname]
    end
    base_un = PBD[varname]
    if length(control_dim.setpts) > 1
        setpts = ustrip.(base_un, control_dim.setpts)
        ramprates = ustrip.(base_un/u"s", control_dim.ramprates)
        holds = ustrip.(u"s", control_dim.holds)
        control_ndim = RampedVariable(setpts, ramprates, holds)
    else
        control_ndim = RampedVariable(ustrip(base_un, control_dim.setpts[1]))
    end
    return control_ndim
end

const PBD = const PARAMS_BASE_DIMS = Dict{Symbol, Any}(
    :Kw => u"W/m^2/K",
    :Kv => u"W/m^2/K",
    :Q_ic => u"W/m^3",
    :Q_gl_RF => u"W",
    :Q_ck => u"W/m^3",
    :k => u"W/m/K",

    :Tsh => u"K",
    :Tw0 => u"K",
    :Tf0 => u"K",
    :Tw => u"K",
    :Tf => u"K",

    :p_ch => u"Pa",

    :kf => u"W/m/K",
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
    :Rp0 => u"m/s",

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
    Kw = 1e2 *u"W/K/m^2" # 1/(Contact resistance 1e-4 + cylindrical resistance from outer wall 1e-3 to 1e-5)
    Kv = 20 * u"W/K/m^2"
    Q_ic = 0.0u"W/cm^3"
    Q_ck = 0.0u"W/m^3"
    m_cp_gl = 5u"g" * LevelSetSublimation.cp_gl # Half of a 10R vial's mass contributing; all of a 2R.

    # Mass transfer
    p_ch = 100u"mTorr" # 100 mTorr is about 13 Pa
    ϵ = 0.9 # 90% porosity
    l = 1e-6u"m" # ~size of a pore
    κ = 1e-10u"m^2" # ~size^2 of a pore
    Rp0 = 1.4u"cm^2*Torr*hr/g" # R0 from Rp: guess from thin-film thickness & pore size?
    k = k_sucrose * (1-ϵ)

    R = 8.3145u"J/mol/K"
    Mw = .018u"kg/mol" #mol/kg
    μ = LevelSetSublimation.μ

    # Sublimation
    ΔH = LevelSetSublimation.ΔH
    # ΔHsub = 678.0 # u"cal/g"
    ρf = ρ_ice 
    Cpf = Cp_ice
    kf = LevelSetSublimation.kf

    # Rw = 8.3145 / .018 # J/molK * mol/kg
    # calc_ρvap(T) = calc_psub(T)/Rw/T


    params = Dict{Symbol, Any}()
    @pack! params = Kw, Kv, Q_ic, Q_ck, k, m_cp_gl, ΔH, kf, ρf, p_ch, ϵ, l, κ, Rp0, R, Mw, μ, Cpf
    # @pack! params = Kw, Kv, k, m_cp_gl, ΔH, kf, ρf, p_ch, ϵ, l, κ, Rp0, R, Mw, μ, Cpf
    return params
end
        
# """
#     make_ramp(ustart, uend, ramprate, ts)

# Convenience function for filling out a time series during and after setpoint ramp.
# `ustart` and `uend` are initial and final setpoints; `ts` the time points at which to sample.

# `ts` need not start at 0, but should contain at least enough time for the full ramp.
# """
# function make_ramp(ustart, uend, ramprate, ts)
#     du = uend - ustart
#     dt = du / ramprate
#     us = fill(ustart, size(ts))
#     dus = min.((ts.-ts[begin])/dt, 1) .* (uend - ustart)
#     return us .+ dus
# end


# --------- Functions mapping state `u` to level set and temperatures

"""
    ϕ_T_from_u_view(u, dom)

Take the current system state `u` and return views corresponding to `ϕ`, `Tf`, and `Tw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u_view(u, dom)
    ϕ = @views reshape(u[1:dom.ntot], size(dom))
    Tf = @view u[dom.ntot+1:dom.ntot+dom.nr]
    Tw = @view u[end]
    return ϕ, Tf, Tw
end

"""
    ϕ_T_from_u(u, dom)

Take the current system state `u` and break it into `ϕ`, `Tf`, and `Tw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u(u, dom)
    ϕ = reshape(u[1:dom.ntot], size(dom))
    Tf = u[dom.ntot+1:dom.ntot+dom.nr]
    Tw = u[end]
    return ϕ, Tf, Tw
end

"""
    make_u0_ndim(init_prof, Tf0, Tw0, dom)

Set up a vector of dynamic variables as initial state for simulation.

Structure of vector `u`:
- `dom.ntot` entries for level set function `ϕ`; initialized with profile using `make_ϕ0`.
- `dom.nr` entries for frozen (ice) temperature `Tf`
- 1 entry for glass (outer wall) temperature `Tw`
"""
function make_u0_ndim(init_prof, Tf0, Tw0, dom)
    Tf0_nd = ustrip.(u"K", Tf0)
    Tw0_nd = ustrip(u"K", Tw0)
    # ϕ0 = make_ϕ0(init_prof, dom)
    # ϕ0_flat = reshape(ϕ0, :)

    ϕ0_flat = reshape(make_ϕ0(init_prof, dom), :)
    u0 = similar(ϕ0_flat, dom.ntot+dom.nr+1) # Add 2 to length: Tf, Tw
    u0[1:dom.ntot] .= ϕ0_flat
    u0[dom.ntot+1:dom.ntot+dom.nr] .= Tf0_nd 
    u0[end] = Tw0_nd
    return u0
end

function make_u0_ndim(config::Dict)
    # dom = Domain(config)
    # Tf0_nd = ustrip.(u"K", config[:Tf0])
    # Tw0_nd = ustrip(u"K", get(config, :Tw0, config[:Tf0]))
    make_u0_ndim(config[:init_prof], config[:Tf0], 
                get(config, :Tw0, config[:Tf0]), Domain(config))
end

# """
#     ϕ_T_into_u!(u, ϕ, Tf, Tw, dom)

# Take `ϕ`, `Tf`, and `Tw`, and stuff them into `u` with appropriate indices.
# Nothing too fancy--just to keep indexing abstract
# """
# function ϕ_T_into_u!(u, ϕ, Tf, Tw, dom)
#     u[1:dom.ntot] = reshape(ϕ, :)
#     u[dom.ntot+dom.nr] = Tf
#     u[end] = Tw
#     return nothing
# end
# """
#     T_into_u!(u, Tf, Tw, dom)

# Take `Tf` and `Tw` and stuff them into `u` with appropriate indices.
# Nothing too fancy--just to keep indexing abstract
# """
# function T_into_u!(u, Tf, Tw, dom)
#     u[dom.ntot+1] = Tf
#     u[dom.ntot+2] = Tw
#     return nothing
# end
