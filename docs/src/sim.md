# Running a Simulation with `LevelSetSublimation`

## Main Simulation Method
```@docs
sim_from_dict
```

## Setting Up for Simulation
```@docs
Domain
make_ϕ0
```

## The Guts

```@docs
LevelSetSublimation.ϕevol_RHS!
LevelSetSublimation.ϕevol_RHS
LevelSetSublimation.reinit_wrap
LevelSetSublimation.next_reinit_time
```