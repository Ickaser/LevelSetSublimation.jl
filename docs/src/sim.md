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

```@doctest
cparams = make_default_params()
init_prof = :circ
Tf0 = 233.15u"K"
QRFvw = 0.002u"W" # = volumetric * relevant vial volume
t_samp = (0:0.1:1) .* u"hr"
Tsh = 263.15u"K"
QRFf = 0.3u"W/cm^3"
pch = 100u"mTorr"

controls = Dict{Symbol, Any}()
@pack! controls = t_samp, QRFvw, Tsh, QRFf, pch


vialsize = "10R"
fillvol = 2u"mL"

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol
```

## The Guts (not exported)

```@docs
dudt_heatmass
dudt_heatmass!
reinit_wrap
```

## Heat-transfer-only simulation

```@docs
sim_heatonly
dudt_heatonly
dudt_heatonly!
```