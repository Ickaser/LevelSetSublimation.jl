
# Plotting and Accessing Simulation Results

## Visualizing simulation results

```@docs
summaryplot
resultsanim
plotframe
```

## Accessing simulation results

The return from `sim_from_dict` is a dictionary, with a field `"dom"` containing a `Domain` and `"sol"` containing an `ODESolution`.

The `ODESolution` type comes from [`DifferentialEquations`](https://docs.sciml.ai/DiffEqDocs/stable/), and comes with lots of convenience functions documented [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/#solution).

```@docs
ϕ_T_from_u
```

At a particular timestep, `u` is used to describe the system's state. 
With the current implementation, `u` is a vector, with `nr*nz` values of the level set function throughout the domain, followed by 1 value for frozen temperature `Tf` and one value for glass temperature `Tvw`.
`ϕ_T_from_u` will return a tuple `(ϕ, Tf, Tvw)` ; this is just so you don't have to hard-code which indices are which variables.


## Fine-grained plotting functions, used internally 

```@docs
plot_cylheat
freshplot
heat
LevelSet
arrows
markfront
```

`plot(::LevelSet, reflect=true)` plots the zero-level set for a level set function `LevelSet(phi, dom::Domain)`. `reflect` reflects across the y-axis to make cylindrical plots.