# Details on the Level Set Method

TODO: fill this out


# Main Level Set Tools

```@docs
reinitialize_ϕ
reinitialize_ϕ!
advect_ϕ
advect_ϕ!
vector_extrap_from_front
```

# Discretizations
```@docs
calc_∇ϕ_1st
calc_Nd∇v
calc_Nd∇v!
𝒢
𝒢_all
calc_Vd∇ϕ
```

# Under-the-Hood Level Set Tools

```@docs
identify_Γ
Γ_cells
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