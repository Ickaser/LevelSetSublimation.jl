# @quickactivate :LevelSetSublimation
const LSS = LevelSetSublimation
using Latexify, TypedTables

# -------------------- Read in LyoPronto data


using CSV
lpfname = datadir("lyopronto", "output_saved_230512_1714.csv")
lpraw = CSV.read(lpfname, Table)
function lyopronto_rename(row)
    nt = (t = row.var"Time [hr]"*u"hr",
        Tsub = row.var"Sublimation Temperature [C]"*u"°C", 
        Tsh = row.var"Shelf Temperature [C]"*u"°C", 
        md = row.var"Sublimation Flux [kg/hr/m^2]"*u"kg/hr/m^2",
        dryfrac = row.var"Percent Dried"/100)
    return nt
end
slice = floor.(Int, range(1, length(lpraw), length=100))
lpdat = map(lyopronto_rename, lpraw[slice])

# Set up simulation

base_props = LSS.base_props

vialsize = "6R"
fillvol = 2u"mL"

simgridsize = (26, 51)
init_prof = :flat 


# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16u"cm*hr*Torr/g"
Tguess = 250u"K"
l = sqrt(u"R"*Tguess/base_props.Mw) / A1
# Other mass transfer params
c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
ϵ = (ρ_solution-c_solid)/ρ_solution 
κ = 0.0u"m^2" # no viscous-regime flow allowed
# Microwave params
B_d = 0.0u"Ω/m^2"
B_f = 0.0u"Ω/m^2"
B_vw = 0.0u"Ω/m^2"
# Kshf value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
ri, ro = LyoPronto.get_vial_radii(vialsize)
Afac = ro^2/ri^2 # Correct for Av/Ap
Kshf = RpFormFit(KC*Afac, KP*Afac, KD)
# Other heat transfer stuff
kd = LSS.k_sucrose * (1-ϵ)
m_v = LyoPronto.get_vial_mass(vialsize)
A_v = π*ro^2
Kvwf = 0.0u"W/m^2/K"
tcprops = TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, B_f, B_vw)

# ---- Variables which can be controlled during a run
f_RF = RampedVariable(0.0u"GHz")
P_per_vial = RampedVariable(0u"W")
Tsh = RampedVariable([238.15u"K", 293.15u"K"], 1u"K/minute")
pch = RampedVariable(150u"mTorr")
tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

paramsd = base_props, tcprops, tvprops

# Pack parameters into a Dict
config = Dict{Symbol, Any}()
@pack! config = paramsd, vialsize, fillvol, simgridsize
config[:time_integ] = Val(:exp_newton)

# Run the simulation
@time res, fname = produce_or_load(sim_from_dict, config; filename=hash, verbose=true, tag=true, prefix=datadir("sims", "lyopronto"));

# @time res = sim_from_dict(config; verbose=true)
# fname = datadir("sims", "lyopronto_"*string(hash(config), base=10)*".jld2")
# safesave(fname, res)

sim = res["sim"];
sim.Tf.u[:,1]
# Get simulation results at matching times
tsol, Tsol, msol, fsol = compare_lyopronto_res(lpdat.t, sim);
# ------------- Comparison plots

# uformatter(l, u) = string(l, raw", $\\left[", latexify(u)[2:end-1], raw"\\right]$")
uf = latexsquareunitlabel

begin
pl1 = plot(u"hr", u"°C", unitformat=uf, ylabel="T_\\text{f}", xlabel="t")
@df lpdat plot!(:t, :Tsub, label="LyoPronto", legend=:bottomright)
plot!(tsol, Tsol, label="LevelSetSublimation")
plot!(Tsh, tmax=6u"hr", label=L"T_\text{sh}", c=:black)
pl2 = plot(u"hr", u"kg/hr/m^2", unitformat=uf, ylabel="\\dot{m}''", xlabel="t")
@df lpdat plot!(:t, :md, label="LyoPronto", unitformat=:square,legend=:bottomright)
plot!(tsol, msol, label="LevelSetSublimation")
pl3 = plot(u"hr", NoUnits, unitformat=uf, ylabel="drying progress", xlabel="t")
@df lpdat plot!(:t, :dryfrac, label="LyoPronto", legend=:bottomright)
plot!(tsol, fsol, label="LevelSetSublimation")
plot(pl1, pl2, pl3)
end

# ------------- Plots for prelim
plot!(pl1, size=(300, 250))
savefig(plotsdir("lyopronto_T.svg"))
savefig(plotsdir("lyopronto_T.pdf"))
plot!(pl2, size=(300, 250))
savefig(plotsdir("lyopronto_subflux.svg"))
savefig(plotsdir("lyopronto_subflux.pdf"))
plot!(pl3, size=(300, 250))
savefig(plotsdir("lyopronto_dryfrac.svg"))
savefig(plotsdir("lyopronto_dryfrac.pdf"))
