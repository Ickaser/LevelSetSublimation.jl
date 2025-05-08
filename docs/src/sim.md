# Running a Simulation with `LevelSetSublimation`

```@meta
CurrentModule = LevelSetSublimation
```

## Main Simulation Method
```@docs
sim_from_dict
```

## Setting Up for Simulation
```@docs
PhysicalProperties
TimeConstantProperties
TimeVaryingProperties
TimeVaryingPropertiesSnapshot
```

Here is a sample simulation setup, in a case where nothing complicated is happening.

```julia
# ---- Properties which do not change in time

# Mass transfer
ϵ = 0.95 # 90% porosity
κ = 0.0u"m^2" # ~size^2 of a pore
Rp0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l = sqrt(base_props.R*Tguess/base_props.Mw) / A1
# Heat transfer
kd = LSS.k_sucrose * (1-ϵ)
m_v = LyoPronto.get_vial_mass(vialsize)
A_v = π*LyoPronto.get_vial_radii(vialsize)[2]^2
# Microwave
B_d = 0.0u"Ω/m^2"
tcprops = TimeConstantProperties(ϵ, l_bulk, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, Bf, Bvw)

# Properties which may change in time
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

paramsd = LSS.base_props, tcprops, tvprops

vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (41, 31)

config = Dict{Symbol, Any}()
@pack! config = paramsd, vialsize, fillvol, simgridsize
```

## Run a simulation

```julia
res = sim_from_dict(config)
```

## Visualize and work with simulation results

## The Guts (not exported)

```@docs
```

## Heat-transfer-only simulation: OUTDATED

```@docs
sim_heatonly
dudt_heatonly
dudt_heatonly!
```