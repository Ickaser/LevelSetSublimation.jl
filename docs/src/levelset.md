# Details on the Level Set Method

TODO: fill this out


# Main Level Set Tools

```@docs
reinitialize_ϕ
reinitialize_ϕ!
reinitialize_ϕ_HCR!
advect_ϕ
advect_ϕ!
extrap_v_pde
extrap_v_fastmarch
```

# Discretizations
```@docs
calc_∇ϕ_1st
calc_Nd∇v
calc_Nd∇v!
𝒢_1st
𝒢_weno
calc_Vd∇ϕ
```

# Under-the-Hood Level Set Tools

```@docs
identify_Γ
identify_B
identify_regions_RC
calc_dϕdr_sdf
calc_dϕdz_sdf
update_ϕ_in_Γ!
```
# Debugging Help
```@docs
plot_RC
compute_frontvel_1
```