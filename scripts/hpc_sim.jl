using Pkg
Pkg.activate("..")
Pkg.instantiate()
using DrWatson
@quickactivate :LevelSetSublimation


# dom = Domain(51, 51, 1.0, 1.0)
cparams = p = make_default_params()

cparams[:Cpf] = 2050u"J/kg/K"
Tf0 = 233.15u"K"

vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (51, 51)

# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
cparams[:l] = l_bulk

init_prof = :flat
Tsh = RampedVariable([233.15u"K", 283.15u"K"], [1u"K/minute"], [10u"hr"])

# With fixed vial mass including whole vial...
cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"
Q_ic = RampedVariable(0.16u"W/cm^3")
Q_gl_RF = RampedVariable(0.13u"W") # = volumetric * relevant vial volume
cparams[:Kw] = 100.0u"W/m^2/K"

p_ch = RampedVariable(100u"mTorr")

controls = Dict{Symbol, Any}()
# @pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol


# # Set up stuff to make debugging easier
params, ncontrols = params_nondim_setup(cparams, controls)

dom = Domain(config)

u0 = make_u0_ndim(init_prof, Tf0, Tf0, dom)
# ϕ0 = make_ϕ0(init_prof, dom)   
# u0 = zeros(dom.ntot+2)
# u0[1:dom.ntot] = reshape(ϕ0, :)
# u0[dom.ntot+1] = ustrip(u"K", Tf0)
# u0[dom.ntot+2] = ustrip(u"K", Tf0)
# reinitialize_ϕ_HCR!(ϕ0, dom)
# T0 = solve_T(u0, dom, params)
# p0 = solve_p(u0, T0, dom, params)

@time res = sim_from_dict(config, verbose=true)

safesave(datadir("sims", "hpcrun.jld2"), res)

resultsanim(res, config, "hpc_test", seconds_length=10)
resultsanim(res, config, "hpc_test", seconds_length=10, heatvar=:p)
