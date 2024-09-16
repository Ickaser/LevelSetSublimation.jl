export make_base_properties, make_M1_properties
export make_ϕ0, make_ramp
export params_nondim_setup, make_u0_ndim
export ϕ_T_from_u, ϕ_T_from_u_view

export TimeConstantProperties, TimeVaryingProperties

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
        ϕ0 = [0.5*(dom.zmax/dom.rmax)*r + z - 0.6dom.zmax + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    elseif ϕtype == :cornercone
        @warn "The more arcane starting shapes are not well-tested."
        ϕ0 = [0.5*(dom.zmax/dom.rmax)*r + z - 0.5dom.zmax + ϵ 
                for r in dom.rgrid, z in dom.zgrid]
    else
        @error "ArgumentError: Invalid ϕ0 kind to make_ϕ0" ϕtype
    end
    return ϕ0
end
struct PhysicalProperties
    cp_gl 
    ρ_gl
    R 
    Mw 
    μ
    ΔH 
    ρf 
    Cpf 
    kf
    ε0 
    εpp_d
    εpp_f 
    εpp_vw 
end


"""
    make_base_properties()

Return a dictionary of base physical properties in SI units.

These should not generally need to change from simulation to simulation, although worth checking periodically.
"""
function make_base_properties()
            
    
    cp_gl = LevelSetSublimation.cp_gl # Half of a 10R vial's mass contributing; all of a 2R.
    ρ_gl = LevelSetSublimation.ρ_gl

    R = 8.3145u"J/mol/K"
    Mw = .018u"kg/mol" #mol/kg
    μ = LevelSetSublimation.μ

    # Sublimation
    ΔH = LevelSetSublimation.ΔH
    ρf = ρ_ice 
    Cpf = Cp_ice
    kf = LevelSetSublimation.kf

    ε0 = LevelSetSublimation.ε0
    εpp_d = 0.0
    εpp_f = LevelSetSublimation.εpp_f
    εpp_vw = LevelSetSublimation.εpp_vw

    props = PhysicalProperties(cp_gl, ρ_gl, R, Mw, μ, ΔH, ρf, Cpf, kf, ε0, εpp_d, εpp_f, εpp_vw)
    return props
end

const base_props = make_base_properties()

"""
"""
struct TimeConstantProperties
    # Mass transfer
    ϵ # 90% porosity
    l # ~size of a pore
    κ # ~size^2 of a pore
    Rp0 # R0 from Rp: guess from thin-film thickness & pore size?
    # Heat transfer
    kd # dry layer thermal conductivity
    Kvwf # vial wall to frozen layer heat transfer coeff
    m_v # vial mass
    A_rad # radiative area for vial wall-shelf heat transfer
    # Microwave
    B_d # dry layer field strength
    B_f # frozen layer field strength coefficient
    B_vw # vial wall field strength coefficient
end

"""
"""
struct TimeVaryingProperties
    f_RF # RF frequency
    P_per_vial # RF power per vial
    Tsh # Shelf temperature
    pch # Chamber pressure
    Kshf # Heat transfer coefficient as function of pressure
end

struct TimeVaryingPropertiesSnapshot
    f_RF # RF frequency
    P_per_vial # RF power per vial
    Tsh # Shelf temperature
    pch # Chamber pressure
    Kshf # Heat transfer coefficient as function of pressure
end

function (tvp::TimeVaryingProperties)(t)
    snap = TimeVaryingPropertiesSnapshot(
        tvp.f_RF(t),
        tvp.P_per_vial(t),
        tvp.Tsh(t),
        tvp.pch(t),
        tvp.Kshf(tvp.pch(t)),
    )
    return snap
end
function (tup::Tuple{PhysicalProperties, TimeConstantProperties, TimeVaryingProperties})(t)
    return (tup[1], tup[2], tup[3](t))
end


function nondim_controlvar(tvp, varname)
    control_dim = getfield(tvp, varname)
    if varname == :Kshf
        control_ndim = p->ustrip(PBD[varname], tvp.Kshf(p*u"Pa"))
        return control_ndim
    end
    if dimension(control_dim(0u"s")) != dimension(PBD[varname])
        @error "Bad units on potentially time-varying variable." varname control_dim PBD[varname]
    end
    base_un = PBD[varname]
    control_ndim = t->ustrip(base_un, control_dim(t*u"s"))
    return control_ndim
end

function nondim_param(tcp, pk)
    p = getfield(tcp, pk)
    if (length(p) > 1) 
        if any([dimension(PBD[pk])] .!= dimension(p))
            @error "Bad dimensions in passed parameter." pk getfield(tcp, pk) PBD[pk]   
        end
    else
        if (dimension(PBD[pk]) != dimension(p))
            @error "Bad dimensions in passed parameter." pk getfield(tcp, pk) PBD[pk]   
        end
    end
    return ustrip.(PBD[pk], p)
end

"""
    params_nondim_setup(params_dim)
    params_nondim_setup(base_d, tcp_d, tvp_d)
    

"""
function params_nondim_setup(base_d, tcp_d, tvp_d)
    nd_b(key) = nondim_param(base_d, key)
    nd_p(key) = nondim_param(tcp_d, key)
    nd_c(key) = nondim_controlvar(tvp_d, key)

    base_n = PhysicalProperties(map(nd_b, fieldnames(typeof(base_d)))...)
    tcp_n = TimeConstantProperties(map(nd_p, fieldnames(typeof(tcp_d)))...)
    tvp_n = TimeVaryingProperties(map(nd_c, fieldnames(typeof(tvp_d)))...)

    return base_n, tcp_n, tvp_n
end
function params_nondim_setup(params_dim)
    base_d, tcp_d, tvp_d = params_dim
    params_nondim_setup(base_d, tcp_d, tvp_d)
end

const PBD = const PARAMS_BASE_DIMS = Dict{Symbol, Any}(
    # Base properties
    :cp_gl => u"J/kg/K",
    :ρ_gl=> u"kg/m^3",
    :R => u"J/mol/K",
    :Mw => u"kg/mol",
    :μ => u"Pa*s",
    :ΔH => u"J/kg",
    :ρf => u"kg/m^3",
    :Cpf => u"J/kg/K",
    :kf => u"W/m/K",
    :ε0 => u"F/m",
    :εpp_d => NoUnits,
    :εpp_f => NoUnits,
    :εpp_vw => NoUnits, 

    # Time constant properties
    :ϵ =>NoUnits, # 90% porosity
    :l =>u"m", # ~size of a pore
    :κ =>u"m^2", # ~size^2 of a pore
    :Rp0 => u"m/s", # R0 from Rp: guess from thin-film thickness & pore size?
    :kd => u"W/m/K",
    :Kvwf => u"W/m^2/K",
    :m_v => u"kg", # vial mass
    :A_rad =>u"m^2", # radiative area for vial wall-shelf heat transfer
    :B_d => u"Ω/m^2", # dry layer field strength
    :B_f =>u"Ω/m^2", # frozen layer field strength coefficient
    :B_vw =>u"Ω/m^2", # vial wall field strength coefficient

    # Time varying properties
    :f_RF =>u"Hz", # RF frequency
    :P_per_vial =>u"W", # RF power per vial
    :Tsh =>u"K", # Shelf temperature
    :pch =>u"Pa", # Chamber pressure
    :Kshf =>u"W/m^2/K", # Heat transfer coefficient as function of pressure

    # Initial conditions
    :Tvw0 => u"K",
    :Tf0 => u"K",
    :Tvw => u"K",
    :Tf => u"K",
)


"""
    make_M1_properties()

Returns dicts of physical parameters for the mannitol experimental case which was used to originally develop and validate this model.

Useful for testing.
"""
function make_M1_properties()
    # ---- Properties which do not change in time
    # Mass transfer
    ϵ = 0.95 # 90% porosity
    l = 1e-6u"m" # ~size of a pore
    κ = 0.0u"m^2" # ~size^2 of a pore
    Rp0 = 1.4u"cm^2*Torr*hr/g" # R0 from Rp: guess from thin-film thickness & pore size?
    # Heat transfer
    kd = k_sucrose * (1-ϵ)
    Kvwf = 24.7 *u"W/K/m^2" 
    m_v = LyoPronto.get_vial_mass("6R")
    A_rad = π*LyoPronto.get_vial_radii("6R")[2]^2
    # Microwave
    B_d = 0.0u"Ω/m^2"
    B_f = 2.8e8u"Ω/m^2"
    B_vw = 4.6e6u"Ω/m^2"

    tcprops = TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_rad, B_f, B_vw)

    # ----- Properties which may change in time
    f_RF = RampedVariable(8.0u"GHz")
    pch = RampedVariable(100.0u"mTorr")
    Tsh = RampedVariable(uconvert.(u"K", [-40.0, 10]*u"°C"), 1u"K/minute")
    P_per_vial = RampedVariable(10u"W"/17)

    # Heat transfer coefficient as function of pressure
    KC = 2.75e-4u"cal/s/K/cm^2"
    KP = 8.93e-4u"cal/s/K/cm^2/Torr"
    KD = 0.46u"1/Torr"
    Kshf = RpFormFit(KC, KP, KD)
    tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

    return base_props, tcprops, tvprops
end
        

# --------- Functions mapping state `u` to level set and temperatures

"""
    ϕ_T_from_u_view(u, dom)

Take the current system state `u` and return views corresponding to `ϕ`, `Tf`, and `Tvw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u_view(u, dom)
    ϕ = @views reshape(u[1:dom.ntot], size(dom))
    Tf = @view u[dom.ntot+1:dom.ntot+dom.nr]
    Tvw = @view u[end]
    return ϕ, Tf, Tvw
end

"""
    ϕ_T_from_u(u, dom)

Take the current system state `u` and break it into `ϕ`, `Tf`, and `Tvw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u(u, dom)
    ϕ = reshape(u[1:dom.ntot], size(dom))
    Tf = u[dom.ntot+1:dom.ntot+dom.nr]
    Tvw = u[end]
    return ϕ, Tf, Tvw
end

"""
    make_u0_ndim(init_prof, Tf0, Tvw0, dom)

Set up a vector of dynamic variables as initial state for simulation.

Structure of vector `u`:
- `dom.ntot` entries for level set function `ϕ`; initialized with profile using `make_ϕ0`.
- `dom.nr` entries for frozen (ice) temperature `Tf`
- 1 entry for glass (outer wall) temperature `Tvw`
"""
function make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    Tf0_nd = ustrip.(u"K", Tf0)
    Tvw0_nd = ustrip(u"K", Tvw0)
    # ϕ0 = make_ϕ0(init_prof, dom)
    # ϕ0_flat = reshape(ϕ0, :)

    ϕ0_flat = reshape(make_ϕ0(init_prof, dom), :)
    u0 = similar(ϕ0_flat, dom.ntot+dom.nr+1) # Add 2 to length: Tf, Tvw
    u0[1:dom.ntot] .= ϕ0_flat
    u0[dom.ntot+1:dom.ntot+dom.nr] .= Tf0_nd 
    u0[end] = Tvw0_nd
    return u0
end

function make_u0_ndim(config::Dict)
    # dom = Domain(config)
    # Tf0_nd = ustrip.(u"K", config[:Tf0])
    # Tvw0_nd = ustrip(u"K", get(config, :Tvw0, config[:Tf0]))
    make_u0_ndim(config[:init_prof], config[:Tf0], 
                get(config, :Tvw0, config[:Tf0]), Domain(config))
end

# """
#     ϕ_T_into_u!(u, ϕ, Tf, Tvw, dom)

# Take `ϕ`, `Tf`, and `Tvw`, and stuff them into `u` with appropriate indices.
# Nothing too fancy--just to keep indexing abstract
# """
# function ϕ_T_into_u!(u, ϕ, Tf, Tvw, dom)
#     u[1:dom.ntot] = reshape(ϕ, :)
#     u[dom.ntot+dom.nr] = Tf
#     u[end] = Tvw
#     return nothing
# end
# """
#     T_into_u!(u, Tf, Tvw, dom)

# Take `Tf` and `Tvw` and stuff them into `u` with appropriate indices.
# Nothing too fancy--just to keep indexing abstract
# """
# function T_into_u!(u, Tf, Tvw, dom)
#     u[dom.ntot+1] = Tf
#     u[dom.ntot+2] = Tvw
#     return nothing
# end
