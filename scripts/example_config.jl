dom = Domain(51, 51, 1.0, 1.0)
params = p = make_artificial_params()
ϕ0type = :circ

config = Dict{Symbol, Any}()
@pack! config = dom, params, ϕ0type

ϕ0 = make_ϕ0(ϕ0type, dom)
reinitialize_ϕ_HCR!(ϕ0, dom)
T0 = solve_T(ϕ0, dom, params)