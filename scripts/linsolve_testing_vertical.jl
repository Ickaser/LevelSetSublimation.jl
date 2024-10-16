const LSS = LevelSetSublimation

using SpecialFunctions
using Roots
using Unitful
using Plots
using FastGaussQuadrature


# charfun_r(l) = besselj0(l*Ri)*bessely1(l*R) - bessely0(l*Ri)*besselj1(l*R)
# λj = [find_zero(charfun_r, (j-0.5)*π/(R-Ri)) for j in 1:400]
# plot(λj)
# ll = range(100, 1000, length=50)
# plot(ll, charfun_r.(ll))
# plot(λj[2:end] .- λj[1:end-1])
# hline!([2.5π/R])


Rp0 = ustrip(u"m/s", 1.4u"Torr*cm^2/g*hr")
A1 = 16u"cm*hr*Torr/g"
Mw = 18u"g/mol"
Rg = 8.3145u"J/mol/K"
T0 = 250.15u"K"
l_bulk = sqrt(Rg*T0/Mw) / A1
b = ustrip(u"s", l_bulk*sqrt(Mw/Rg/T0)) #Mw/R/T * (l*NaNMath.sqrt(R*T/Mw)
Δp = ustrip(u"Pa", LSS.calc_psub(T0) - 100u"mTorr")



# -------------------

cparams = make_default_params()
cparams[:l] = upreferred(l_bulk)
cparams[:κ] = 0.0u"m^2"
cparams[:Rp0] = Rp0*u"m/s"

cparams[:Kshf] *= 0
cparams[:Kvwf] = 10.0u"W/m^2/K"

simgridsize = (101, 101)
Tfm = fill(ustrip(u"K", T0), simgridsize[1])
Tdm = fill(ustrip(u"K", T0), simgridsize)
Tf0 = T0
Tvw0 = T0 + 20u"K"

controls = Dict{Symbol, Any}()
# @pack! controls = t_samp, QRFvw, Tsh, QRFf, pch
QRFvw = RampedVariable(0.0u"W")
Tsh = RampedVariable(0.0u"K")
QRFf = RampedVariable(0.0u"W/cm^3")
pch = RampedVariable(100u"mTorr")
@pack! controls = QRFvw, Tsh, QRFf, pch

init_prof = :flat
vialsize = "6R"
fillvol = 5u"mL"
config = Dict{Symbol, Any}()
@pack! config = simgridsize, cparams, init_prof, Tf0, Tvw0, controls, vialsize, fillvol

dom = Domain(config)
params, ncontrols = params_nondim_setup(cparams, controls)
um = LSS.make_u0_ndim(config)
ϕm = @views reshape(um[iϕ(dom)], size(dom))
ϕm .+= .8*dom.zmax - 1e-6 + 1e-8

R = dom.rmax
L = dom.zmax
hl = LSS.get_subf_z(ϕm, dom)

# LSS.eval_b(Tdm, Tdm, params)



p_num = solve_p(um, Tfm, Tdm, dom, params)
# p_sol1 = gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=250)
p_sol(z, h) = -Δp/(b*Rp0 + L-h)*(z-h)
p_anl = [LSS.calc_psub(ustrip(u"K", T0)) + (z < hl ? 0 : p_sol(z, hl))
            for r in dom.rgrid, z in dom.zgrid]
heat(p_num, dom)
pl1 = heat(abs.(p_num .- p_anl)./p_anl, dom)

function gen_vertprof(er)
    um = LSS.make_u0_ndim(config)
    ϕm = @views reshape(um[iϕ(dom)], size(dom))
    ϕm .+= .8dom.zmax + er - 0.99e-6#-2e-6

    R = dom.rmax
    L = dom.zmax
    hl = LSS.get_subf_z(ϕm, dom)

    p_num = solve_p(um, Tfm, Tdm, dom, params)
    p_anl = [LSS.calc_psub(ustrip(u"K", T0)) + (z < hl ? 0 : p_sol(z, hl))
                for r in dom.rgrid, z in dom.zgrid]
    return [p_num[1,:], p_anl[1,:]]
end

perturbs = range(0, 2dom.dz, length=100) 
perturbs = perturbs .- step(perturbs)
@time res = [gen_vertprof(ei) for ei in perturbs]
plot(res[1][1])
plot([ri[1][22] for ri in res])
plot!([ri[2][22] for ri in res])
plot(perturbs./dom.dz, [abs(ri[1][42]-ri[2][42])/ri[2][42] for ri in res])
plot(hcat(res...)')
res

T_num = solve_T(um, Tfm, dom, params)



plot(dom.rgrid, p_num[:,end], label="numerical")
plot!(dom.rgrid, p_anl[:,end], label="analytical", c=:red)
vline!([Ri])
plot(p_num[:,end] .- p_anl[:,end])
vline!([11, 12])
plot!(xlim=(Ri*0.8, Ri*1.2))

plot( p_num[:,end-1], label="numerical")
plot!(p_anl[:,end-1], label="analytical", c=:red)
plot(p_num[:,end-1] .- p_anl[:,end-1])
# # Set up stuff to make debugging easier
# rr = range(0.005, .01, length=100)
# zz = range(0, L, length=110)
# @time ps = [p_sol1(r, z) for z in zz, r in rr]
# heatmap(rr, zz, ps, cmap=:plasma)

