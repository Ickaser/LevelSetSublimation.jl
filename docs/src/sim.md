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
make_default_params
```

Here is a sample simulation setup, in a case where nothing complicated is happening.

```@doctest
cparams = make_default_params()
ϕ0type = :circ
Tf0 = 233.15u"K"
Q_gl_RF = 0.002u"W" # = volumetric * relevant vial volume
t_samp = (0:0.1:1) .* u"hr"
Tsh = 263.15u"K"
Q_ic = 0.3u"W/cm^3"
p_ch = 100u"mTorr"

controls = Dict{Symbol, Any}()
@pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch


vialsize = "10R"
fillvol = 2u"mL"

config = Dict{Symbol, Any}()
@pack! config = cparams, ϕ0type, Tf0, controls, vialsize, fillvol
```

## The Guts (not exported)

```@docs
uevol_heatmass
uevol_heatmass!
reinit_wrap
next_reinit_time
```

## Heat-transfer-only simulation

```@docs
sim_heatonly
uevol_heatonly
uevol_heatonly!
```