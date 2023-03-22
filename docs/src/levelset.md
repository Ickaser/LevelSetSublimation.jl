# Details on the Level Set Method

TODO: fill this out


## Main Level Set Tools

```@docs
reinitialize_ϕ_HCR!
reinitialize_ϕ_HCR
extrap_v_fastmarch
```

## Discretizations, Backend

```@docs
𝒢_weno  
𝒢_weno_all  
identify_Γ
Γ_cells
wenodiffs_local
```

## Obsolete Level Set Tools

```@docs
reinitialize_ϕ
reinitialize_ϕ!
advect_ϕ
advect_ϕ!
extrap_v_pde
```

## Obsolete Discretizations, etc.

```@docs
calc_Vd∇ϕ
calc_Nd∇v
calc_Nd∇v!
𝒢_1st  
calc_∇ϕ_1st
compute_frontvel_1
plot_RC
identify_regions_RC
calc_dϕdr_sdf
calc_dϕdz_sdf
update_ϕ_in_Γ!
```


## Under-the-Hood Level Set Tools

```@docs
identify_B
```
## Debugging Help
```@docs
```