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
Q_gl_RF = RampedVariable(0.08u"W") # = volumetric * relevant vial volume
cparams[:Kw] = 100.0u"W/m^2/K"

p_ch = RampedVariable(100u"mTorr")

controls = Dict{Symbol, Any}()
# @pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol


@time res = sim_from_dict(config, verbose=true)


@time Tf_sol = virtual_thermocouple([0, 0, 0.1, 1], [0.1, 0, 0, 0], res, config)

resconf = @strdict res config Tf_sol

safesave(datadir("sims", "hpcrun.jld2"), resconf)

# resultsanim(res, config, "hpc_test", seconds_length=10)
# resultsanim(res, config, "hpc_test", seconds_length=10, heatvar=:p)
