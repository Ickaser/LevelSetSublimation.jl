export make_M1_properties
export PhysicalProperties
export TimeConstantProperties
export TimeVaryingProperties
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
        ϕ0 = [z - dom.zmax*0.5+ϵ for r in dom.rgrid, z in dom.zgrid]
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


"""
An immutable type for representing a slate of physical properties and constants needed.

A default constructor is provided which uses all defaults (borosilicate glass, water, ice), but thanks to [`Parameters.jl`](https://mauro3.github.io/Parameters.jl/stable/)
you can call with just values that need to be changed, e.g.

    PhysicalProperties() # all defaults
    PhysicalProperties(ρ_vw = 1u"kg/m^3", cp_vw = 1u"J/kg/K", εpp_vw = 1e-1) # adjust vial wall properties

All properties stored here are:
$(FIELDS)
"""
@with_kw struct PhysicalProperties
    "Universal gas constant."
    R = 8.3145u"J/mol/K"
    "Vacuum permittivity (universal constant)."
    ε0 = LevelSetSublimation.ε0
    "Density of vial wall (defaults to borosilicate glass)."
    ρ_vw = LevelSetSublimation.ρ_gl
    "Heat capacity of vial wall (defaults to borosilicate glass)."
    cp_vw = LevelSetSublimation.cp_gl # Half of a 10R vial's mass contributing; all of a 2R.
    "Dielectric loss coefficient of vial wall (defaults to borosilicate glass)."
    εpp_vw = LevelSetSublimation.εpp_gl
    "Density of frozen material (taken as water ice)."
    ρf = ρ_ice 
    "Heat capacity of frozen material (taken as water ice)."
    Cpf = cp_ice
    "Thermal conductivity of frozen material."
    kf = LevelSetSublimation.kf
    "Molecular weight of sublimating species (defaults to water)."
    Mw = .018u"kg/mol" 
    "Gas phase viscosity (in rarefied regime) of sublimating species (defaults to water)."
    μ = LevelSetSublimation.μ
    "Heat of sublimation of sublimating species (defaults to water); give as positive number."
    ΔH = LevelSetSublimation.ΔH
    "Dielectric loss coefficient of frozen layer (defaults to water ice)."
    εppf = LyoPronto.εppf
    "Dielectric loss coefficient of dry layer (defaults to 0)."
    εpp_d = 0.0
end


"""
An instance of [`PhysicalProperties`](@ref) with default values.
"""
const base_props = PhysicalProperties()

"""
    TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, B_f, B_vw)

A struct for holding physical properties which are likely to change from case to case.

No default constructor is provided by intention--all of these parameters should be at least considered before running a simulation.

$(FIELDS)
"""
struct TimeConstantProperties
    # Mass transfer
    "dimensionless, the porosity or fraction of solid space taken up by ice."
    ϵ 
    "length, a dusty gas model parameter (roughly the pore size)"
    l 
    "area, a dusty gas model parameter (roughly the area of a pore neck)"
    κ 
    "length/time, empirical mass transfer resistance at zero dry layer height"
    Rp0 
    # Heat transfer
    "dry layer thermal conductivity"
    kd 
    "vial wall to frozen layer heat transfer coeff"
    Kvwf 
    "vial mass"
    m_v 
    "total vial-bottom area, used for heat transfer"
    A_v 
    # Microwave
    "Ω/m^2, dry layer field strength coefficient"
    B_d 
    "Ω/m^2, frozen layer field strength coefficient"
    B_f 
    "Ω/m^2, vial wall field strength coefficient"
    B_vw 
end

"""
    TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

A struct for holding controlled parameters that may change over time.

To get the value of all those parameters at a given time, call the struct 
with a time, and it will evaluate each field at that time and provide a [`TimeVaryingPropertiesSnapshot`](@ref). 

No default constructor is provided by intention--all of these parameters should be at least considered before running a simulation.

Each should be passed as a callable, which returns the value of the parameter as a function of time, with the exception of `Kshf` which is a function of pressure.
See the `RampedVariable` and `RpFormFit` types from `LyoPronto`, which are intended to make this more convenient.

$(FIELDS)
"""
struct TimeVaryingProperties
    "Microwave frequency"
    f_RF 
    "Microwave power per vial"
    P_per_vial 
    "Shelf temperature"
    Tsh 
    "Chamber pressure"
    pch 
    "Heat transfer coefficient as function of pressure"
    Kshf 
end

"""
    TimeVaryingPropertiesSnapshot(f_RF, P_per_vial, Tsh, pch, Kshf)

A struct for holding a snapshot of the controlled parameters that may change over time.

This is meant to be constructed by calling an instance of the [`TimeVaryingProperties`](@ref) type.

$(FIELDS)
"""
struct TimeVaryingPropertiesSnapshot
    "Microwave frequency"
    f_RF 
    "Microwave power per vial"
    P_per_vial 
    "Shelf temperature"
    Tsh 
    "Chamber pressure"
    pch 
    "Heat transfer coefficient as function of pressure"
    Kshf 
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

const SimSetup = Tuple{PhysicalProperties, TimeConstantProperties, TimeVaryingProperties} 
function (tup::SimSetup)(t)
    return (tup[1], tup[2], tup[3](t))
end


struct NondimensionalizedFunc{F, I, O}
    f::F
    in_un::I
    out_un::O
end
(nf::NondimensionalizedFunc)(x...) = ustrip(nf.out_un, nf.f((x .* nf.in_un)...))

LyoPronto.extract_ts(nd::NondimensionalizedFunc) = LyoPronto.extract_ts(nd.f, un=nd.in_un)
function LyoPronto.get_tstops(tvp::TimeVaryingProperties)
    return get_tstops((tvp.f_RF, tvp.P_per_vial, tvp.Tsh, tvp.pch))
end
LyoPronto.get_tstops(tup::SimSetup) = LyoPronto.get_tstops(tup[3])

function nondim_controlvar(tvp, varname)
    control_dim = getfield(tvp, varname)
    if varname == :Kshf
        # control_ndim = p->ustrip(PBD[varname], tvp.Kshf(p*u"Pa"))
        # return control_ndim
        return NondimensionalizedFunc(tvp.Kshf, u"Pa", PBD[varname])
    end
    if dimension(control_dim(0u"s")) != dimension(PBD[varname])
        @error "Bad units on potentially time-varying variable." varname control_dim PBD[varname]
    end
    base_un = PBD[varname]
    # control_ndim = t->ustrip(base_un, control_dim(t*u"s"))
    return NondimensionalizedFunc(control_dim, u"s", base_un)
end

function nondim_param(tcp, pk)
    p = getfield(tcp, pk)
    if pk == :εppf
        var_ndim = (T,f)->ustrip(PBD[pk], tcp.εppf(T*PBD[:Tsh], f*PBD[:f_RF]))
        var_ndim = NondimensionalizedFunc(tcp.εppf, (PBD[:Tsh], PBD[:f_RF]), PBD[pk])
        return var_ndim
    end
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
    :cp_vw => u"J/kg/K",
    :ρ_vw=> u"kg/m^3",
    :R => u"J/mol/K",
    :Mw => u"kg/mol",
    :μ => u"Pa*s",
    :ΔH => u"J/kg",
    :ρf => u"kg/m^3",
    :Cpf => u"J/kg/K",
    :kf => u"W/m/K",
    :ε0 => u"F/m",
    :εpp_d => NoUnits,
    :εpp_vw => NoUnits, 
    :εppf => NoUnits,

    # Time constant properties
    :ϵ =>NoUnits, # 90% porosity
    :l =>u"m", # ~size of a pore
    :κ =>u"m^2", # ~size^2 of a pore
    :Rp0 => u"m/s", # R0 from Rp: guess from thin-film thickness & pore size?
    :kd => u"W/m/K",
    :Kvwf => u"W/m^2/K",
    :m_v => u"kg", # vial mass
    :A_v =>u"m^2", # radiative area for vial wall-shelf heat transfer
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

Returns physical parameters for the mannitol experimental case which was used to originally develop and validate this model.

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
    A_v = π*LyoPronto.get_vial_radii("6R")[2]^2
    # Microwave
    B_d = 0.0u"Ω/m^2"
    B_f = 2.8e8u"Ω/m^2"
    B_vw = 4.6e6u"Ω/m^2"

    tcprops = TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, B_f, B_vw)

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

# The source of truth for what indices are used for what variables!!!
"""
    iϕ(dom::Domain) = 1:dom.ntot

This is the source of truth for what indices in the system state are used for `ϕ`.
See also [`iTf`](@ref) and [`iTvw`](@ref).
"""
iϕ(dom::Domain) = 1:dom.ntot
"""
    iTf(dom::Domain) = dom.ntot+1:dom.not+nr

This is the source of truth for what indices in the system state are used for `Tf`.
See also [`iϕ`](@ref) and [`iTvw`](@ref).
"""
iTf(dom::Domain) = dom.ntot+1:dom.ntot+dom.nr
"""
    iTvw(dom::Domain) = dom.ntot+dom.nr+1 # returns a single index
    iTvw(dom::Domain, dummyarg) = dom.ntot+dom.nr+1:dom.ntot+dom.nr+1 # returns a 1-length UnitRange

This is the source of truth for what indices in the system state are used for `Tvw`.
See also [`iϕ`](@ref) and [`iTf`](@ref).
"""
iTvw(dom::Domain) = dom.ntot+dom.nr+1
iTvw(dom::Domain, dummyarg) = dom.ntot+dom.nr+1:dom.ntot+dom.nr+1
"""
    ulen(dom::Domain) = iTvw(dom)

Returns the total length of the system state vector, which currently is the same as the index of `Tvw`.
"""
ulen(dom::Domain) = iTvw(dom)

ϕ_from_u(u, dom) = reshape(u[iϕ(dom)], size(dom))
# reshape(a::AbstractArray, dom::Domain) = reshape(a, size(dom))

"""
    ϕ_from_u_view(u, dom)

Take the current system state `u` and return views corresponding to `ϕ`, `Tf`, and `Tvw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u_view(u, dom)
    ϕ = @views reshape(u[iϕ(dom)], size(dom))
    Tf = @view u[iTf(dom)]
    Tvw = @view u[iTvw(dom)]
    return ϕ, Tf, Tvw
end

"""
    ϕ_T_from_u(u, dom)

Take the current system state `u` and break it into `ϕ`, `Tf`, and `Tvw`.
Nothing too fancy--just to avoid rewriting the same logic everywhere
"""
function ϕ_T_from_u(u, dom)
    ϕ = reshape(u[iϕ(dom)], size(dom))
    Tf = u[iTf(dom)]
    Tvw = u[iTvw(dom)]
    return ϕ, Tf, Tvw
end

"""
    make_u0_ndim(init_prof, Tf0, Tvw0, dom)

Set up a vector of dynamic variables as initial state for simulation.

Structure of vector `u`: set by [`iϕ`](@ref), [`iTf`](@ref), and [`iTvw`](@ref).
"""
function make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    Tf0_nd = ustrip.(u"K", Tf0)
    Tvw0_nd = ustrip(u"K", Tvw0)

    ϕ0_flat = reshape(make_ϕ0(init_prof, dom), :)
    u0 = similar(ϕ0_flat, ulen(dom)) 
    u0[iϕ(dom)] .= ϕ0_flat
    u0[iTf(dom)] .= Tf0_nd 
    u0[iTvw(dom)] = Tvw0_nd
    return u0
end

function make_u0_ndim(config::Dict)
    # dom = Domain(config)
    Tf0_nd = ustrip.(u"K", get(config,:Tf0, config[:paramsd].Tsh(0u"s")))
    Tvw0_nd = ustrip(u"K", get(config, :Tvw0, Tf0_nd*u"K"))
    make_u0_ndim(config[:init_prof], Tf0_nd, Tvw0_nd, Domain(config))
end