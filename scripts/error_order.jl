@quickactivate :LevelSetSublimation
const LSS = LevelSetSublimation
using UnitfulLatexify

# plot defaults
default(:linewidth, 3)
default(:markersize, 5)
default(:fontfamily, "Computer Modern")

# -----------------------
# Set up physical parameters
begin

base_props = LSS.base_props

# Get some stuff from LyoProntoNIIMBLRF
lcfitparams = load(datadir("exp_pro", "M1_KvRpRF.jld2"))
@unpack Kvwf, Bf, Bvw = lcfitparams
expdat = load(datadir("exp_pro", "M1_processed.jld2"))
thmdat = expdat["thm_pd"]
t_end = expdat["t_end"]

vialsize = "6R"
fillvol = 5u"mL"
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

f_RF = RampedVariable(8.0u"GHz")
pch = RampedVariable(100.0u"mTorr")
Tsh = RampedVariable(uconvert.(u"K", [-40.0, 10]*u"°C"), 1u"K/minute")
P_per_vial = RampedVariable(10u"W"/17 * 0.54)
# Heat transfer coefficient as function of pressure
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

paramsd = LSS.base_props, tcprops, tvprops
end

# ------------
# Assemble simulation
simgridsizes = [(21, 15),
                (31, 21),
                (31, 31),
                (41, 31),
                (51, 31),
                (61, 31),
                (71, 31),
                (31, 41),
                (31, 51),
                (31, 61),
                (31, 71),
                (41, 41),
                (51, 51),
                (61, 61),
                (71, 71),
]

# config = Dict{Symbol, Any}()
# @pack! config = paramsd, vialsize, fillvol, simgridsize
config = @dict paramsd vialsize fillvol
config[:simgridsize] = simgridsizes
config[:time_integ] = :dae_then_exp
# config[:time_integ] = :exp_newton
allconfigs = dict_list(config)


# Run simulations
timing = map(allconfigs) do config
    @elapsed produce_or_load(sim_from_dict, config; filename=hash, verbose=true, tag=false, prefix=datadir("sims", "err"))
end

#15-element Vector{Float64}:
#   206.9773286
#    48.5550164
#    65.7792836
#   179.6836682
#   269.2749256
#   342.0172003
#   916.4128808
#   142.8181078
#   217.4506176
#   246.0629959
#   340.6295969
#   225.4222626
#   562.6863109
#  1651.0493648
#  3891.9363282

using DataFrames
allsims = collect_results(datadir("sims"), rinclude=[r"err"])
simtab = Table(map(eachrow(allsims)) do row
    # if row.sim.config[:time_integ] == :dae_then_exp
    #     row.sim = merge(row.sim, (Tf = nothing,))
    # end
    merge(row.sim, (path=row["path"],), NamedTuple(row.sim.config))
end)

best = simtab[argmax(prod.(simtab.simgridsize))];
rbest = simtab[findmax(x->x[1], simtab.simgridsize)[2]];
zbest = simtab[findmax(x->x[2], simtab.simgridsize)[2]];
best.simgridsize

locs = [(0.0, 0.0)]#, (0.8, 0.05), (0.2, 0.5)]
t = range(0, best.sol.t[end]*0.99, length=100)
Tbest = virtual_thermocouple(locs, t, best)
Tbestr = virtual_thermocouple(locs, t, rbest)
Tbestz = virtual_thermocouple(locs, t, rbest)
function compare_T(sim; Tb = Tbest, t=t)
    Ts = virtual_thermocouple(locs, t, sim)
    return sqrt(sum(abs2, Ts-Tb))
end

sqerr = map(simtab) do sim
    lerr = @time compare_T(sim)
    rerr = compare_T(sim, Tb=Tbestr)
    zerr = compare_T(sim, Tb=Tbestz)
    terr = abs(sim.sol.t[end] - best.sol.t[end])
    return (sqerr=lerr, rerr=rerr, zerr=zerr, terr=terr, tend=sim.sol.t[end], nr=sim.simgridsize[1], nz=sim.simgridsize[2])
end
# @time compare_T(simtab[1])
errtab = Table(simtab, sqerr);
sort!(errtab, by=x->x.nr*x.nz);

cs = cgrad(:viridis)
begin
blankplothrC(ylabel="Temperature difference to reference")
# plot!(t, Tbest[:,1], c=cs[end], label="Reference")
markers=[:circle, :heptagon, :utriangle, :dtriangle, :ltriangle, :rtriangle]
for (i,sim) in enumerate(errtab)
    c = cs[257 - sim.nr*sim.nz*256÷prod(best.simgridsize)]
    mark = markers[(i-1)%6+1]
    plot!(t, Tbest[:,1] .- virtual_thermocouple([locs[1]], t, sim), c=c, marker=mark, label="nr=$(sim.nr), nz=$(sim.nz)")
    scatter!([sim.tend], [0], label="", c=c)
end
plot!(xlim=(0e4, 5e4), ylim=(-0.2,0.2), legend=:outerright)
end
plot!(ylim=(-1.9, 0.3))

resetfontsizes()
scalefontsizes(1.2)
begin
plr = @df filter(x->(x.nz==31 && x.sqerr>0), errtab) scatter(:nr, :sqerr, scale=:log10, marker=:square, label=L"n_z=31")
@df filter(x->(x.nr==x.nz && x.sqerr>0), errtab) scatter!(:nr, :sqerr, scale=:log10, label=L"n_z=n_r")
# @df filter(x->(x.nr==31 && x.sqerr>0), errtab) scatter!(:nr, :sqerr, scale=:log10, label="nr=31")
plot!([10^1.5, 10^1.8], [10^0.25, 10^-0.35], label="slope -2", ls=:dash, c=:gray)
plot!([10^1.4, 10^1.9], [10^0.27, 10^0.27], label="slope 0", ls=:dot, c=:gray)
plot!(xlabel=L"$n_r$, number of $r$ grid points", ylabel=L"L2 error in $T_\mathrm{f}$ [K]")

plz = @df filter(x->(x.nr==31 && x.sqerr>0), errtab) scatter(:nz, :sqerr, scale=:log10, marker=:square, label=L"n_r=31")
@df filter(x->(x.nr==x.nz && x.sqerr>0), errtab) scatter!(:nz, :sqerr, scale=:log10, label=L"n_r=n_z")
plot!([10^1.3, 10^1.9], [10^0.8, 10^-0.4], label="slope -2", ls=:dash, c=:gray)
plot!(xlabel=L"$n_z$, number of $z$ grid points", )#ylabel=L"L2 error in $T_\mathrm{f}$ [K]")

plot!(plr, left_margin=20Plots.px, )
plot(plr, plz, layout=(1,2), size=(800, 400))
plot!(legend=:bottomleft, bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_convergence_dae.svg"))
savefig(plotsdir("M1_LSS_convergence_dae.pdf"))

begin
plr = @df filter(x->(x.nz==31 && x.sqerr>0), errtab) scatter(:nr, :terr, scale=:log10, marker=:square, label=L"n_z=31")
@df filter(x->(x.nr==x.nz && x.sqerr>0), errtab) scatter!(:nr, :terr, scale=:log10, label=L"n_z=n_r")
# @df filter(x->(x.nr==31 && x.sqerr>0), errtab) scatter!(:nr, :sqerr, scale=:log10, label="nr=31")
plot!([10^1.5, 10^1.8], [10^2.24, 10^1.34], label="slope -3", ls=:dash, c=:gray)
# plot!([10^1.4, 10^1.9], [10^1.27, 10^1.27], label="slope 0", ls=:dot, c=:gray)
plot!(xlabel=L"$n_r$, number of $r$ grid points", ylabel=L"L2 error in $t_\mathrm{end}$ [s]")

plz = @df filter(x->(x.nr==31 && x.sqerr>0), errtab) scatter(:nz, :terr, scale=:log10, marker=:square, label=L"n_r=31")
@df filter(x->(x.nr==x.nz && x.sqerr>0), errtab) scatter!(:nz, :terr, scale=:log10, label=L"n_r=n_z")
plot!([10^1.35, 10^1.85], [10^2.27, 10^2.27], label="slope 0", ls=:dot, c=:gray)
plot!([10^1.5, 10^1.8], [10^2.24, 10^1.34], label="slope -3", ls=:dash, c=:gray)
plot!(xlabel=L"$n_z$, number of $z$ grid points", )#ylabel=L"L2 error in $T_\mathrm{f}$ [K]")

plot!(plr, left_margin=20Plots.px, )
plot(plr, plz, layout=(1,2), size=(800, 400))
plot!(legend=:bottomleft, bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_convergence_dae_tend.svg"))
savefig(plotsdir("M1_LSS_convergence_dae_tend.pdf"))