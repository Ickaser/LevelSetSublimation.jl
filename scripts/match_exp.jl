using CSV

fname_sim1 = "2023-10-24.jld2"

sim1 = load(datadir("sims", fname_sim1))
@unpack res, config, Tf_sol = sim1
@unpack sol, dom = res

simt = range(0, sol.t[end], length=100)
Tf_sol

plot(simt[1:end-1]./3600, Tf_sol[1:end-1,:] .-273.15, labels=permutedims(["TC$i" for i in 1:5]))
plot!(sol.t[1:end-5]./3600, sol[end,1:end-5] .-273.15, label="Tw")

plot!(xlabel="time (hr)", ylabel="temperature (°C)")
savefig(plotsdir("2023-10-24_multiT.svg"))


dat1 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "TemperatureLog_030223.csv"))
dat2 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "Mannitol_5ml_Temp.csv"), comment="#", stripwhitespace=true)
start_ti = findfirst(dat2["SPLYO.VACUUM_SP.F_CV"] .== 100)
# Pch = Pch[start_ti:end]
Pch = dat2["SPLYO.CHAMBER_CM.F_CV"][start_ti:end]
Tsh = dat2["SPLYO.SHELF_SP.F_CV"][start_ti:end]
Pch_pir = dat2["SPLYO.CHAMBER_PIRANI.F_CV"][start_ti:end]
plT = plot(Tsh)
plp = plot(Pch, ylim=(0, 300))
plp = plot!(Pch_pir)
plot(plT, plp, layout=(2,1))

# start_ti = 92
# endti = 4400

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

plot( tdat, T4dat, label="v1, T4")
plot!(tdat, T1dat, label="v2, T1")
plot!(tdat, T3dat, label="gl, T3")
plot!(tdat, Tsh_rv.(tdat), label="shelf")


# for fname in [fname_sim1, fname_sim2, fname_sim3]
# stuff = map([fname_sim2, fname_sim3]) do fname
#     sim = load(datadir("sims", fname))
#     t = sim["res"]["sol"].t
#     # T = @time virtual_thermocouple(sim["res"], sim["config"])
#     # (sim, t, T)
# end

