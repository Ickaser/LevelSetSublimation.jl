# Details on the Level Set Method

TODO: fill this out

```@meta
CurrentModule = LevelSetSublimation
```

## Main Level Set Tools

```@docs
reinitialize_ϕ_HCR!
reinitialize_ϕ_HCR
extrap_v_fastmarch!
```

## Discretizations, Backend

```@docs
identify_Γ
Γ_cells
𝒢_weno  
𝒢_weno_all  
dϕdx_all_WENO
wenodiffs_local
```

## Error evaluation

```@docs
identify_B
sdf_err_L1
sdf_err_L∞
```

## Volume and surface integrals and evaluations

```@docs
compute_icevol_H
compute_discrete_H
compute_local_H
compute_icesurf_δ
compute_discrete_δ
compute_local_δ
```

## Initial sublimation front shapes, for numerical testing

```@docs
make_ϕ0
```
