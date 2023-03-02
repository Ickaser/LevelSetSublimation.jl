dom = Domain(25, 21, 1.0, 1.0)
T_params = p = make_artificial_params()
ϕ0type = :circ
sim_dt = 1.0

config = Dict{Symbol, Any}()
@pack! config = dom, T_params, ϕ0type, sim_dt

ϕ0 = make_ϕ0(ϕ0type, dom)
reinitialize_ϕ!(ϕ0, dom)
T0 = solve_T(ϕ0, dom, T_params)