@quickactivate :LevelSetSublimation
const LSS = LevelSetSublimation
using UnitfulLatexify

# plot defaults
default(:linewidth, 3)
default(:markersize, 5)

# -----------------------

base_props = LSS.base_props

# Get some stuff from LyoProntoNIIMBLRF
lcfitparams = load(datadir("exp_pro", "M3_KvRpRF.jld2"))
@unpack Kvwf, Bf, Bvw = lcfitparams
expdat = load(datadir("exp_pro", "M3_processed.jld2"))
P_per_vial = expdat["P_per_vial"]
thmdat = expdat["thm_pd"]
t_end = expdat["t_end"]

# ---- Properties which do not change in time
vialsize = "6R"
fillvol = 3u"mL"
# Mass transfer
c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
ϵ = (ρ_solution-c_solid)/ρ_solution 
κ = 0.0u"m^2" # no viscous-regime flow allowed
# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l = sqrt(base_props.R*Tguess/base_props.Mw) / A1
# Heat transfer
kd = LSS.k_sucrose * (1-ϵ)
m_v = LyoPronto.get_vial_mass(vialsize)
A_v = π*LyoPronto.get_vial_radii(vialsize)[2]^2
# Microwave
B_d = 0.0u"Ω/m^2"
tcprops = TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, Bf, Bvw)

# ------- Properties which may change in time
f_RF = RampedVariable(18.0u"GHz")
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable(258.15u"K")
# P_per_vial loaded with data above

# Heat transfer coefficient as function of pressure
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

# paramsd = LSS.base_props, tcprops, tvprops
paramsd = LSS.base_props, tcprops, tvprops

# ------------
# Assemble simulation
simgridsize = (61, 61)

# config = Dict{Symbol, Any}()
# @pack! config = paramsd, vialsize, fillvol, simgridsize
config = @dict paramsd vialsize fillvol simgridsize
config[:time_integ] = :dae_then_exp

# Run simulation

@time res, fname = produce_or_load(sim_from_dict, config; filename=hash, verbose=true, tag=true, prefix=datadir("sims", "M3"))

# @time res = sim_from_dict(config; verbose=true)
# fname = datadir("sims", "M1_"*string(hash(config), base=10)*".jld2")
# safesave(fname, res)

# fname = datadir("sims", "M1_"*string(hash(config), base=10)*".jld2")
# res = load(fname)["sim"]

sim = res["sim"]
@unpack dom = sim

# --------------

# Plot things of interest

# ---------------------------
# Plots

# summaryT(sim, tstart=0.01, tend=0.99, layout=(2,3), clims=(-30, 50))
# plot!(size=(800,300), bottom_margin=10Plots.px)
# savefig(plotsdir("M1_LSS_summary.svg"))
# savefig(plotsdir("M1_LSS_summary.pdf"))

locs = [(0.0, 0.0), (0.8, 0.05), (0.2, 0.5)]
vtmarks = [:diamond, :utriangle, :heptagon]

resetfontsizes()
scalefontsizes(1.2)
begin
pl_sum = summaryT(sim, tstart=0.01, tend=0.80, layout=(1,3))
placethermocouples!(dom, locs, c=palette(:Oranges_4)[4:-1:2], markers=vtmarks, label="", markersize=8);
plot!(size=(800, 200), left_margin=-5Plots.px, right_margin=0Plots.px)
plot!(pl_sum[4], cbar_title="\nTemperature [°C]", right_margin=20Plots.px)
end
savefig(plotsdir("M3_LSS_summVT.svg"))
savefig(plotsdir("M3_LSS_summVT.pdf"))

resetfontsizes()
begin
labs = vcat([L"$T$, model "*i for i in ["bottom","corner","center"]], L"$T_\mathrm{vw}$, model")
pl_T = blankplothrC(;xlabel=L"t", )
plot!(Tsh, c=:black, tmax=t_end, label=L"T_\mathrm{sh}")
@df thmdat exptfplot!(:t, :T1, :T4, lw=3)
@df thmdat exptvwplot!(:t, :T3, lw=3, thickness_scaling=1.0)
vt_plot!(sim, locs; labels=permutedims(labs), markers=permutedims(vtmarks), step=40)
tendplot!(t_end)
plot!(legend=:outerright, size=(800, 400), bottom_margin=15Plots.px, left_margin=15Plots.px)
end
savefig(plotsdir("M3_multiT.svg"))
savefig(plotsdir("M3_multiT.pdf"))
# plot!(pl_sum, margin_right=26Plots.px)
plot(pl_T, pl_sum, layout=@layout([a{0.5h}; b]), size=(500,400))
savefig(plotsdir("M3_allLSS.svg"))
savefig(plotsdir("M3_allLSS.pdf"))

# savefig(plotsdir("M1_multiT_small.svg"))
# savefig(plotsdir("M1_multiT_small.pdf"))


begin
pl, T = plotframe(10*60*60, sim)
placethermocouples!(sim.dom, locs, c=palette(:Oranges_4)[4:-1:2], markers=vtmarks, label="", markersize=8);
# texts = [text(L"$T_\textrm{f%$i}$", 10) for i in 1:3]
# labelthermocouples!(sim.dom, locs, texts)
plot!(xlabel="", cbar=nothing)
# annotate!(pl, dom.rmax.*locx, dom.zmax.*locy, texts)
plot!(size=(240,200), colorbar_title="T [°C]")
end
# savefig(pl, plotsdir("M1_virtual_thermocouple_small.svg"))
# savefig(pl, plotsdir("M1_virtual_thermocouple_small.pdf"))

# --------- Compute effective Rp

t_Rp, hd_eff, Rp_eff = get_eff_Rp(sim);
Rp_orig = @. Rp0 + A1*hd_eff
relerr = (Rp_eff .- Rp_orig)./Rp_orig
begin
pl1 = plot(u"cm", u"cm^2*Torr*hr/g", ylabel="R_p", xlabel="h_d", unitformat=latexsquareunitlabel)
plot!(hd_eff, Rp_eff, label="LSS with curvature")
plot!(hd_eff, Rp_orig, label="Original Rp")
pl2 = plot(u"hr", NoUnits, ylabel="Relative Difference", xlabel="Time")
plot!(t_Rp, relerr; ylim=(-0.1, 0.1))
plot(pl1, pl2, layout=(2,1))
end
savefig(plotsdir("M3_curvatureRp.svg"))
savefig(plotsdir("M3_curvatureRp.pdf"))

# -------- LCLyo model

A_p, A_v = π.* get_vial_radii("6R").^2
c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
LCparams = [
    (RpFormFit(Rp0, A1, 0u"1/cm"), fillvol/A_p, c_solid, ρ_solution),
    (Kshf, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (fillvol*ρ_solution, base_props.Cpf, m_v, base_props.cp_vw),
    (f_RF(0), base_props.εpp_f, base_props.εpp_vw),
    (Kvwf, Bf, Bvw)
]
LCParams = ParamObjRF(LCparams)

lcprob = ODEProblem(LCParams)
lclyo = solve(lcprob, Rodas3())

begin
blankplothrC()
@df thmdat exptfplot!(:t, :T1, :T4, marker=:circle, )
@df thmdat exptvwplot!(:t, :T3, marker=:square, )
modrftplot!(lclyo, labels=permutedims(["\$T_f\$, model "*i for i in ["bottom","side","top"]]))
plot!(size=(400, 200), legend=:outerright)
end
savefig(plotsdir("M3_lc_small.svg"))
savefig(plotsdir("M3_lc_small.pdf"))