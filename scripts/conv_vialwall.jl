# @quickactivate :LevelSetSublimation

# In this script, we allow heat exchange with the vial wall--otherwise, this behaves like LyoPronto default

cparams = make_default_params()

Tf0 = -35.857u"°C"
vialsize = "6R"
fillvol = 2u"mL"

simgridsize = (51, 51)
init_prof = :flat 


# ---- Variables which can be controlled during a run

Tsh = RampedVariable([238.15u"K", 293.15u"K"], [1u"K/minute"], [10u"hr"])
p_ch = RampedVariable(150u"mTorr")
Q_gl_RF = RampedVariable(0u"W")
Q_ic = RampedVariable(0u"W/cm^3")

# Kv value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
cparams[:Kv] = (KC + KP*p_ch(0u"s")/(1+KD*p_ch(0u"s"))) 
cparams[:Kv] *= 3.8/3.14 # Correct for Av/Ap factor

# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0
A1 = 16u"cm*hr*Torr/g"
Tguess = 250u"K"
l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
l_bulk = upreferred(l_bulk)
cparams[:l] = l_bulk
cparams[:ϵ] = 0.95 

cparams[:κ] *= 0 # Multiply by 0, to match dimensions


# Assemble parameters into objects -------------------------

controls = Dict{Symbol, Any}()
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol

# config[:dudt_func] = LSS.dudt_heatmass_dae!

# -------------------- Read in LyoPronto data


# -------------------- Run simulation
@time res = sim_from_dict(config)

# ------------- Comparison plots
# pl1 = plot(tsol, [T_lp, Tsol], ylabel="Tp", labels=permutedims(["LyoPronto", "LevelSetSublimation"]), legend=:bottomright)
# pl2 = plot(tsol, [m_lp, msol], ylabel="sub. flux", labels=permutedims(["LyoPronto", "LevelSetSublimation"]), legend=:bottomright)
# pl3 = plot(tsol, [dryfrac_lp, fsol], ylabel="drying progress", labels=permutedims(["LyoPronto", "LevelSetSublimation"]))

# tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)
# pl1 = plot(t_lp, T_lp, label="LyoPronto", ylabel="Tp", legend=:bottomright)
# plot!(tsol, Tsol, label="LevelSetSublimation")
# pl2 = plot(t_lp, m_lp, label="LyoPronto", ylabel="sub. flux", legend=:bottomright)
# plot!(tsol, msol, label="LevelSetSublimation")
# pl3 = plot(t_lp, dryfrac_lp, label="LyoPronto", ylabel="drying progress", legend=:bottomright)
# plot!(tsol, fsol, label="LevelSetSublimation")
# plot(pl1, pl2, pl3)