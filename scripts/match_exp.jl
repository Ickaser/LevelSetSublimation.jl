using CSV

# fname_sim1 = "2023-10-24.jld2"

# sim1 = load(datadir("sims", fname_sim1))
# @unpack res, config, Tf_sol = sim1
# @unpack sol, dom = res

# simt = range(0, sol.t[end], length=100)
# Tf_sol

# plot(simt[1:end-1]./3600, Tf_sol[1:end-1,:] .-273.15, labels=permutedims(["TC$i" for i in 1:5]))
# plot!(sol.t[1:end-5]./3600, sol[end,1:end-5] .-273.15, label="Tw")

# plot!(xlabel="time (hr)", ylabel="temperature (°C)")
# savefig(plotsdir("2023-10-25_multiT.svg"))


dat1 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "TemperatureLog_030223.csv"))
dat2 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "Mannitol_5ml_Temp.csv"), comment="#", stripwhitespace=true)
start_ti = findfirst(dat2["SPLYO.VACUUM_SP.F_CV"] .== 100)
# Pch = Pch[start_ti:end]
Pch = dat2["SPLYO.CHAMBER_CM.F_CV"][start_ti:end]
Tsh = dat2["SPLYO.SHELF_SP.F_CV"][start_ti:end]
Pch_pir = dat2["SPLYO.CHAMBER_PIRANI.F_CV"][start_ti:end]

# plT = plot(Tsh)
# plp = plot(Pch, ylim=(0, 300))
# plp = plot!(Pch_pir)
# plot(plT, plp, layout=(2,1))


tdat_raw = dat1["Elapsed [s]"]*u"s"
T1dat_raw = dat1["T1 [C]"]*u"°C"
T3dat_raw = dat1["T3 [C]"]*u"°C"
T4dat_raw = dat1["T4 [C]"]*u"°C"
sti = 92
# fti = length(tdat_raw)
fti = argmax(T4dat_raw)
tdat = tdat_raw[sti:fti]
T1dat = T1dat_raw[sti:fti]
T3dat = T3dat_raw[sti:fti]
T4dat = T4dat_raw[sti:fti]


tdat .-= tdat[1]
tdat = uconvert.(u"hr", tdat)


@show tdat[end] - tdat[begin]

Tsh_rv = RampedVariable(uconvert.(u"K", [-40u"°C", 10u"°C"]).+0.0u"K", [0.5u"K/minute"], [1u"hr"])



sim = load(datadir("sims", "2024-04-11.jld2"))
sim = load(datadir("sims", "M1_10658494824944108369.jld2"))

sol = sim["res"]["sol"]
# @unpack Tf_sol = sim
# Tf_sol *= u"K"
# plot!(tsol, Tf_sol, c=2)
tsol = range(0, sol.t[end], 50) *u"s"
locx, locy = [0.0, 0.8, 0.1, 0.5], [0.00, 0.1, 0.5, 0.9]
Tf_sol = virtual_thermocouple(locx, locy, ustrip.(u"s", tsol), sim["res"], sim["config"]).*u"K"
Tvw_m = sol.(ustrip.(u"s", tsol), idxs=2653).*u"K"

default(:fontfamily, "Helvetica")
default(:palette, :tab20c)
resetfontsizes()
scalefontsizes(1.5)
trim = range(1, length(tdat), step=20)
begin
plot( tdat, T4dat, label=L"T_\textrm{f1}"*", exp.", c=1)
plot!(tdat, T1dat, label=L"T_\textrm{f2}"*", exp.", c=2)
plot!(tdat[trim], T3dat[trim], label=L"T_\textrm{vw}"*", exp.", c=1, ls=:dash)
plot!(tdat, Tsh_rv.(tdat), label="shelf", c=:black)
labs = permutedims(["\$T_\\textrm{f$i}\$, model" for i in 1:4])
plot!(tsol, Tf_sol, c=permutedims(5:8), label=labs)
plot!(tsol, Tvw_m, c=5, ls=:dash, label=L"T_\textrm{vw}"*", model",)
plot!(xlabel="Time [hr]", ylabel="Temperature [°C]")
plot!(legend_font="sans-serif", legend_column=2)
vline!([380/36], c=:gray, linealpha=0.5, lw=5, label="")
vline!([12, 12.1], c=:gray, ls=[:dash, :dot], label="")
plot!(size=(500,400))
end
savefig(plotsdir("multiT.svg"))

pl, T = plotframe(38000.0, sim["res"], sim["config"])
pl;
texts = [text("\$T_\\textrm{f$i}\$", "Computer Modern") for i in 1:4]
dom = sim["res"]["dom"]
plot()
plot!(xlabel="");
scatter!(pl, dom.rmax.*locx, dom.zmax.*locy, c=5:8, label="", markersize=20);
annotate!(pl, dom.rmax.*locx, dom.zmax.*locy, texts);
savefig(pl, plotsdir("virtual_thermocouple.svg"))

# for fname in [fname_sim1, fname_sim2, fname_sim3]
# stuff = map([fname_sim2, fname_sim3]) do fname
#     sim = load(datadir("sims", fname))
#     t = sim["res"]["sol"].t
#     # T = @time virtual_thermocouple(sim["res"], sim["config"])
#     # (sim, t, T)
# end

