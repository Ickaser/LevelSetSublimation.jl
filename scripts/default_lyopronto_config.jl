# @quickactivate :LevelSetSublimation

cparams = make_default_params()

Tf0 = -35.857u"°C"
vialsize = "6R"
fillvol = 2u"mL"

simgridsize = (51, 51)
init_prof = :flat 


# ---- Variables which can be controlled during a run
# Set points

t_samp = range(0, 2, step=1//60) .* u"hr"
# Ramp from -35 to 20
Tsh = fill(20.0u"°C", length(t_samp))
Tsh[1:101] .= make_ramp(-35u"°C", 20u"°C", 1u"K/minute", t_samp[1:101])
p_ch = cparams[:p_ch] = 150u"mTorr"
Q_gl_RF = 0u"W"
Q_ic = 0u"W/cm^3"

# Kv value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
cparams[:Kv] = (KC + KP*p_ch/(1+KD*p_ch)) 
cparams[:Kv] *= 3.8/3.14 # Correct for Av/Ap factor

# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0
A1 = 16u"cm*hr*Torr/g"
Tguess = 240u"K"
l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
l_bulk = upreferred(l_bulk)
cparams[:l] = l_bulk
cparams[:ϵ] = 0.95 

cparams[:κ] *= 0 # Multiply by 0, to match dimensions
cparams[:Kgl] *= 0

# Assemble parameters into objects -------------------------

controls = Dict{Symbol, Any}()
@pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol


# -------------------- Read in LyoPronto data

using CSV
lpfname = datadir("lyopronto", "output_saved_230512_1714.csv")
lpdat = CSV.File(lpfname)

slice = floor.(Int, range(1, length(lpdat), length=100))
t_lp = lpdat["Time [hr]"][slice]*u"hr"
T_lp = lpdat["Sublimation Temperature [C]"][slice]*u"°C"
Tsh_lp = lpdat["Shelf Temperature [C]"][slice]*u"°C"
m_lp = lpdat["Sublimation Flux [kg/hr/m^2]"][slice]*u"kg/hr/m^2"
dryfrac_lp = lpdat["Percent Dried"][slice] / 100

# -------------------- Run simulation
@time res = sim_from_dict(config)

# ------------- Comparison plots
tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)
pl1 = plot(tsol, [T_lp, Tsol], ylabel="Tp", labels=permutedims(["LyoPronto", "LevelSetSublimation"]), legend=:bottomright)
pl2 = plot(tsol, [m_lp, msol], ylabel="sub. flux", labels=permutedims(["LyoPronto", "LevelSetSublimation"]))
pl3 = plot(tsol, [dryfrac_lp, fsol], ylabel="drying progress", labels=permutedims(["LyoPronto", "LevelSetSublimation"]))
plot(pl1, pl2, pl3)