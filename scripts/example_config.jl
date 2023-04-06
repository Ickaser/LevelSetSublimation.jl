dom = Domain(51, 51, 1.0, 1.0)
cparams = p = make_artificial_params()
ϕ0type = :box
Tf0 = 233.15

config = Dict{Symbol, Any}()
@pack! config = dom, cparams, ϕ0type, Tf0

ϕ0 = make_ϕ0(ϕ0type, dom)   
u0 = zeros(dom.ntot+1)
u0[1:dom.ntot] = reshape(ϕ0, :)
u0[dom.ntot+1] = Tf0
reinitialize_ϕ_HCR!(ϕ0, dom)
T0 = solve_T(u0, dom, cparams)
p0 = solve_p(u0, T0, dom, cparams)