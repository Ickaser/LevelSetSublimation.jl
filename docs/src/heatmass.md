# Solving steady-state heat equation


```@meta
CurrentModule = LevelSetSublimation
```

TODO: fill this out

## Physical Equation Solution
```@docs
solve_T
solve_p
solve_p_given_b
```

## Lumped computations

```@docs
compute_Qice
compute_topmassflux
compute_Qice_noflow
compute_Qice_nodry
```


## Computing velocity: coupling solutions to level set motion

```@docs
compute_frontvel_mass
compute_frontvel_heat
compute_frontvel_fixedspeed
```

## Geometric computations from level set

```@docs
compute_icesh_area
compute_icegl_area
```

See also [`compute_icevol_H`](@ref) and [`compute_icesurf_δ`](@ref).