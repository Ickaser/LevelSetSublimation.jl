# @quickactivate :LevelSetSublimation

cparams = make_default_params()

Tf0 = -35.857u"°C"
vialsize = "6R"
fillvol = 2u"mL"

simgridsize = (26, 51)
init_prof = :flat 


# ---- Variables which can be controlled during a run

Tsh = RampedVariable([238.15u"K", 293.15u"K"], [1u"K/minute"], [])
p_ch = RampedVariable(150u"mTorr")
QRFvw = RampedVariable(0u"W")
QRFf = RampedVariable(0u"W/cm^3")

# Kshf value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
cparams[:Kshf] = (KC + KP*p_ch(0u"s")/(1+KD*p_ch(0u"s"))) 
cparams[:Kshf] *= 3.8/3.14 # Correct for Av/Ap factor

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
cparams[:Kvwf] *= 0

# Assemble parameters into objects -------------------------

controls = Dict{Symbol, Any}()
@pack! controls = QRFvw, Tsh, QRFf, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol

# config[:dudt_func] = LSS.dudt_heatmass_dae!

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
@time sim, fname = produce_or_load(sim_and_postprocess, config; filename=hash, verbose=false, tag=true, prefix=datadir("sims", "lyopronto"))
# @time sim = sim_and_postprocess(config)
# safesave(datadir("sims", savename(fi)))

res = sim["res"]

# ------------- Comparison plots
# pl1 = plot(tsol, [T_lp, Tsol], ylabel="Tp", labels=permutedims(["LyoPronto", "LevelSetSublimation"]), legend=:bottomright)
# pl2 = plot(tsol, [m_lp, msol], ylabel="sub. flux", labels=permutedims(["LyoPronto", "LevelSetSublimation"]), legend=:bottomright)
# pl3 = plot(tsol, [dryfrac_lp, fsol], ylabel="drying progress", labels=permutedims(["LyoPronto", "LevelSetSublimation"]))

tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)
pl1 = plot(t_lp, T_lp, label="LyoPronto", ylabel=L"T_\textrm{f}", unitformat=:square, legend=:bottomright)
plot!(tsol, Tsol, label="LevelSetSublimation")
pl2 = plot(t_lp, m_lp, label="LyoPronto", legend=:bottomright)
plot!(tsol, msol, label="LevelSetSublimation")
plot!(ylabel=P"$\dot{m}''$  [$kg\ hr^{-1} m^{-2}$]", )
pl3 = plot(t_lp, dryfrac_lp, label="LyoPronto", ylabel="drying progress", legend=:bottomright)
plot!(tsol, fsol, label="LevelSetSublimation")
plot(pl1, pl2, pl3)

# ------------- Plots for prelim
plot!(pl1, size=(300, 250))
savefig(plotsdir("lyopronto_T.pdf"))
plot!(pl2, size=(300, 250))
savefig(plotsdir("lyopronto_subflux.pdf"))
plot!(pl3, size=(300, 250))
savefig(plotsdir("lyopronto_dryfrac.pdf"))
