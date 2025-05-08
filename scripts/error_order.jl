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
l = sqrt(u"R"*Tguess/base_props.Mw) / A1
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
simgridsizes = [
                (21, 21),
                (31, 31),
                (41, 41),
                (51, 51),
                (61, 61),
                (71, 71),
                # (41, 31),
                # (51, 31),
                # (61, 31),
                # (71, 31),
                # (31, 21),
                # (31, 41),
                # (31, 51),
                # (31, 61),
                # (31, 71),
]

# config = Dict{Symbol, Any}()
# @pack! config = paramsd, vialsize, fillvol, simgridsize
config = @dict paramsd vialsize fillvol
config[:simgridsize] = simgridsizes
config[:time_integ] = Val(:dae_then_exp)
# config[:time_integ] = Val(:exp_newton)
allconfigs = dict_list(config)


# Run simulations
timing = map(allconfigs) do config
    @elapsed produce_or_load(sim_from_dict, config; filename=hash, verbose=true, tag=false, prefix=datadir("sims", "err"))
end

#15-element Vector{Float64}: # out of order now
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
    # if row.sim.config[:time_integ] == Val(:dae_then_exp)
    #     row.sim = merge(row.sim, (Tf = nothing,))
    # end
    merge(row.sim, (path=row["path"],), NamedTuple(row.sim.config))
end)

best = simtab[argmax(prod.(simtab.simgridsize))];
rbest = simtab[findmax(x->x[1], simtab.simgridsize)[2]];
zbest = simtab[findmax(x->x[2], simtab.simgridsize)[2]];
best.simgridsize

locs = [(0.0, 0.0)]#, (0.8, 0.05), (0.2, 0.5)]
t1 = range(0, 40000, length=80) # stop before ice detaches from wall
t2 = range(42000, best.sol.t[end]*0.99, length=20) # pick up after
Tbest1 = virtual_thermocouple(locs, t1, best)
Tbest2 = virtual_thermocouple(locs, t2, best)
# Tbestr = virtual_thermocouple(locs, t, rbest)
# Tbestz = virtual_thermocouple(locs, t, rbest)
function compare_Tf(sim, Tb, t)
    Ts = virtual_thermocouple(locs, t, sim)
    return sqrt(sum(abs2, Ts-Tb)/length(Ts))
end
function compare_Tvw(sim, best, t)
    Ts = sim.sol.(t).Tvw
    Tb = best.sol.(t).Tvw
    return sqrt(sum(abs2, Ts-Tb)/length(Ts))
end

simtab_rz = filter(x->x.simgridsize[1]==x.simgridsize[2], simtab)

errs = map(simtab_rz) do sim
    err_Tf_rz1 = @time compare_Tf(sim, Tbest1, t1)
    err_Tf_rz2 = @time compare_Tf(sim, Tbest2, t2)
    # err_Tf_r = compare_Tf(sim, Tb=Tbestr)
    # err_Tf_z = compare_Tf(sim, Tb=Tbestz)
    err_t_rz = abs(sim.sol.t[end] - best.sol.t[end])/best.sol.t[end]
    err_td_rz = abs(sim.sol.tsplit - best.sol.tsplit)/best.sol.tsplit
    # err_t_r = abs(sim.sol.t[end] - rbest.sol.t[end])
    # err_t_z = abs(sim.sol.t[end] - zbest.sol.t[end])
    err_Tvw_rz1 = compare_Tvw(sim, best, t1)
    err_Tvw_rz2 = compare_Tvw(sim, best, t2)
    return merge((nr=sim.simgridsize[1], nz=sim.simgridsize[2]), @ntuple err_Tf_rz1 err_Tf_rz2 err_Tvw_rz1 err_Tvw_rz2 err_t_rz err_td_rz)
end
# @time compare_Tf(simtab[1])
errtab = Table(simtab_rz, errs);
sort!(errtab, by=x->x.nr*x.nz);

cs = cgrad(:viridis)
begin
pl1 = blankplothrC(ylabel="Temperature difference to reference")
# plot!(t, Tbest[:,1], c=cs[end], label="Reference")
markers=[:circle, :heptagon, :utriangle, :dtriangle, :ltriangle, :rtriangle]
for (i,sim) in enumerate(errtab)
    c = cs[257 - sim.nr*sim.nz*256÷prod(best.simgridsize)]
    mark = markers[(i-1)%6+1]
    plot!(t1/3600, Tbest1[:,1] .- virtual_thermocouple([locs[1]], t1, sim), c=c, marker=mark, label=L"n_r=n_z=%$(sim.nr)")
end
plot!(xlim=(0, 11.3), ylim=(-0.2,0.2), legend=nothing, left_margin=20Plots.px)
pl2 = blankplothrC(ylabel="")
for (i,sim) in enumerate(errtab)
    c = cs[257 - sim.nr*sim.nz*256÷prod(best.simgridsize)]
    mark = markers[(i-1)%6+1]
    plot!(t2/3600, Tbest2[:,1] .- virtual_thermocouple([locs[1]], t2, sim), c=c, marker=mark, label=L"n_r=n_z=%$(sim.nr)")
end
plot!(legend=:outerright)
plot(pl1, pl2, size=(800, 400), bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_error_over_time.svg"))
savefig(plotsdir("M1_LSS_error_over_time.pdf"))

resetfontsizes()
scalefontsizes(1.2)
begin
# plr = @df filter(x->(x.nz==31 && x.rerr>0), errtab) scatter(:nr, :rerr, scale=:log10, marker=:square, label=L"n_z=31")
# @df filter(x->(x.nr==x.nz && x.rzerr>0), errtab) scatter!(:nr, :rzerr, scale=:log10, label=L"n_z=n_r")
# # @df filter(x->(x.nr==31 && x.rzerr>0), errtab) scatter!(:nr, :rzerr, scale=:log10, label="nr=31")
# plot!([10^1.5, 10^1.8], [10^0.25, 10^-0.35], label="slope -2", ls=:dash, c=:gray)
# plot!([10^1.4, 10^1.9], [10^0.01, 10^0.01], label="slope 0", ls=:dot, c=:gray)
# plot!(xlabel=L"$n_r$, number of $r$ grid points", ylabel=L"L2 error in $T_\mathrm{f}$ [K]")

# plz = @df filter(x->(x.nr==31 && x.zerr>0), errtab) scatter(:nz, :zerr, scale=:log10, marker=:square, label=L"n_r=31")
# @df filter(x->(x.nr==x.nz && x.rzerr>0), errtab) scatter!(:nz, :rzerr, scale=:log10, label=L"n_r=n_z")
# plot!([10^1.3, 10^1.9], [10^0.8, 10^-0.4], label="slope -2", ls=:dash, c=:gray)
# plot!([10^1.4, 10^1.9], [10^0.31, 10^0.31], label="slope 0", ls=:dot, c=:gray)
# plot!(xlabel=L"$n_z$, number of $z$ grid points", )#ylabel=L"L2 error in $T_\mathrm{f}$ [K]")

plTf = @df filter(x->(x.err_Tf_rz1>0), errtab) scatter(:nr, :err_Tf_rz1, scale=:log10, label="pre-detachment")
@df filter(x->(x.err_Tf_rz1>0), errtab) scatter!(:nr, :err_Tf_rz2, scale=:log10, marker=:square, label="post-detachment")
# @df filter(x->(x.nr==31 && x.rzerr>0), errtab) scatter!(:nr, :rzerr, scale=:log10, label="nr=31")
plot!([10^1.3, 10^1.8], [10^0.2, 10^-0.8], label=L"err$\propto n^{-2}$", ls=:dash, c=:gray)
plot!([10^1.3, 10^1.8], [10^-0.9, 10^-1.9], label="", ls=:dash, c=:gray)
# plot!([10^1.4, 10^1.9], [10^0.01, 10^0.01], label="slope 0", ls=:dot, c=:gray)
plot!(xlabel=L"$n_r=n_z$, number of grid points each direction", ylabel=L"RMS error in $T_\mathrm{f}$ [K]")

plot!(legend=:bottomleft, bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_convergence_dae_Tf.svg"))
savefig(plotsdir("M1_LSS_convergence_dae_Tf.pdf"))

begin
plTvw = @df filter(x->(x.err_Tvw_rz1>0), errtab) scatter(:nr, :err_Tvw_rz1, scale=:log10, label="pre-detachment")
@df filter(x->(x.err_Tvw_rz1>0), errtab) scatter!(:nr, :err_Tvw_rz2, scale=:log10, marker=:square, label="post-detachment")
# @df filter(x->(x.nr==31 && x.rzerr>0), errtab) scatter!(:nr, :rzerr, scale=:log10, label="nr=31")
plot!([10^1.3, 10^1.8], [10^0.22, 10^-0.78], label=L"err$\propto n^{-2}$", ls=:dash, c=:gray)
plot!([10^1.3, 10^1.8], [10^-0.55, 10^-1.55], label="", ls=:dash, c=:gray)
# plot!([10^1.3, 10^1.8], [10^-0.9, 10^-1.4], label=L"err$\propto n^{-1}$", ls=:dot, c=:gray)
plot!(xlabel=L"$n_r=n_z$, number of grid points each direction", ylabel=L"RMS error in $T_\mathrm{vw}$ [K]")
plot!(legend=:bottomleft, bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_convergence_dae_Tvw.svg"))
savefig(plotsdir("M1_LSS_convergence_dae_Tvw.pdf"))

begin
# plr = @df filter(x->(x.nz==31 && x.rerr>0), errtab) scatter(:nr, :trerr, scale=:log10, marker=:square, label=L"n_z=31")
# @df filter(x->(x.nr==x.nz && x.rzerr>0), errtab) scatter!(:nr, :trzerr, scale=:log10, label=L"n_z=n_r")
# # @df filter(x->(x.nr==31 && x.rzerr>0), errtab) scatter!(:nr, :rzerr, scale=:log10, label="nr=31")
# plot!([10^1.5, 10^1.8], [10^2.24, 10^1.34], label="slope -3", ls=:dash, c=:gray)
# # plot!([10^1.4, 10^1.9], [10^1.27, 10^1.27], label="slope 0", ls=:dot, c=:gray)
# plot!(xlabel=L"$n_r$, number of $r$ grid points", ylabel=L"L2 error in $t_\mathrm{end}$ [s]")
# plot!( legend=:bottomleft)

# plz = @df filter(x->(x.nr==31 && x.tzerr>0), errtab) scatter(:nz, :tzerr, scale=:log10, marker=:square, label=L"n_r=31")
# @df filter(x->(x.nr==x.nz && x.trzerr>0), errtab) scatter!(:nz, :trzerr, scale=:log10, label=L"n_r=n_z")
# plot!([10^1.35, 10^1.85], [10^1.77, 10^1.77], label="slope 0", ls=:dot, c=:gray)
# plot!([10^1.5, 10^1.8], [10^2.24, 10^1.34], label="slope -3", ls=:dash, c=:gray)
# plot!(xlabel=L"$n_z$, number of $z$ grid points", )#ylabel=L"L2 error in $T_\mathrm{f}$ [K]")
# plot!(legend=:bottom)
plt = @df filter(x->x.err_t_rz>0, errtab) scatter(:nz, :err_t_rz, scale=:log10, label=L"t_\mathrm{end}")
@df filter(x->x.err_t_rz>0, errtab) scatter!(:nz, :err_td_rz, scale=:log10, marker=:square, label=L"t_\mathrm{detach}")
# plot!([10^1.35, 10^1.85], [10^1.77, 10^1.77], label="slope 0", ls=:dot, c=:gray)
plot!([10^1.3, 10^1.8], [10^-1.65, 10^-3.15], label=L"err$\propto n^{-3}$", ls=:dash, c=:gray)
plot!(xlabel=L"$n_z=n_r$, number of grid points each direction", ylabel=L"relative error in $t$")
plot!(legend=:bottom)

# plot!(plr, left_margin=20Plots.px, )
# plot(plr, plz, layout=(1,2), size=(800, 400))
plot!(bottom_margin=20Plots.px)
end
savefig(plotsdir("M1_LSS_convergence_dae_t.svg"))
savefig(plotsdir("M1_LSS_convergence_dae_t.pdf"))