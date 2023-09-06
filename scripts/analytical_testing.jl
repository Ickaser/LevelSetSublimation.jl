const LSS = LevelSetSublimation

using SpecialFunctions
using Roots
using Unitful
using Plots
using FastGaussQuadrature


charfun_r(l) = besselj0(l*Ri)*bessely1(l*R) - bessely0(l*Ri)*besselj1(l*R)
λj = [find_zero(charfun_r, (j-0.5)*π/(R-Ri)) for j in 1:400]
plot(λj)
ll = range(100, 1000, length=50)
plot(ll, charfun_r.(ll))
plot(λj[2:end] .- λj[1:end-1])
hline!([2.5π/R])


function gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=100)
    charfun_r(l) = besselj0(l*Ri)*bessely1(l*R) - bessely0(l*Ri)*besselj1(l*R)
    λj = [find_zero(charfun_r, (j-0.5)*π/(R-Ri)) for j in 1:Nmax]
    # λj = zeros(Nmax)
    # λj[1] = find_zero(charfun_r, π/R)
    # for i in 2:Nmax
    #     prev_root = λj[i-1]
    #     @show prev_root
    #     root1 = find_zero(charfun_r, prev_root + 2.5π/R)
    #     if isapprox(root1, prev_root+2.5π/R, atol=10)
    #         λj[i] = root1
    #     else
    #         @warn "problems"
    #     end
    # end
    charfun_z(m) = cos(m*L) - Rp0*b*m*sin(m*L)
    # μk = permutedims([find_zero(charfun_z, (k-0.5)*π/L) for k in 1:Nmax])
    μk = zeros(1, Nmax)
    μk[1] = find_zero(charfun_z, 0.5*π/L)
    for i in 2:Nmax
        prev_root = μk[i-1]
        root1 = find_zero(charfun_z, prev_root + π/L)
        if isapprox(root1, prev_root+π/L, atol=10)
            μk[i] = root1
        else
            @warn "problems"
        end
    end

    function make_uj(j)
        uj(r) = -bessely0(λj[j]*Ri)*besselj0(λj[j]*r) + besselj0(λj[j]*Ri)*bessely0(λj[j]*r)
    end
    function make_duj(j)
        duj(r) = λj[j]* bessely0(λj[j]*Ri)*besselj1(λj[j]*r) - λj[j]*besselj0(λj[j]*Ri)*bessely1(λj[j]*r)
    end
    u = [make_uj(j) for j in 1:Nmax]
    function make_vk(k)
        vk(z) = cos(μk[k]*z)
    end
    v = [make_vk(k) for k in 1:Nmax]

    x, w = gausslegendre(Nmax)
    r_nodes = @. (x +1)/2*(R-Ri) + Ri
    r_wts = @. w/2*(R-Ri)
    z_nodes = @. (x +1)/2*L
    z_wts = @. w/2*L
    innerprod_r(f, g) = sum(@. r_wts * r_nodes * f(r_nodes) * g(r_nodes))
    innerprod_z(f, g) = sum(@. z_wts * f(z_nodes) * g(z_nodes))
    u_normsq = [innerprod_r(uj, uj) for uj in u]
    v_normsq = [innerprod_z(vk, vk) for vk in v]

    # Using p' = p-psub, z-transform then r-transform
    Ck = [vk(L) for vk in v] .* (-Δp/Rp0/b)
    Cjk = [innerprod_r(x->Ci, uj) for uj in u, Ci in Ck] 
    coeffs = Cjk ./(λj.^2 .+ μk.^2) ./ u_normsq ./ permutedims(v_normsq)

    function p_sol(r, z)
        sum(coeffs .* [uj(r) for uj in u] .* permutedims([vk(z) for vk in v]))
    end
    return p_sol
end

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

simgridsize = (101, 101)
Tfm = fill(ustrip(u"K", T0), simgridsize[1])
Tdm = fill(ustrip(u"K", T0), simgridsize)
Tf0 = T0

controls = Dict{Symbol, Any}()
# @pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch
Q_gl_RF = RampedVariable(0.0u"W")
Tsh = RampedVariable(0.0u"K")
Q_ic = RampedVariable(0.0u"W/cm^3")
p_ch = RampedVariable(100u"mTorr")
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch

init_prof = :rad
vialsize = "6R"
fillvol = 5u"mL"
config = Dict{Symbol, Any}()
@pack! config = simgridsize, cparams, init_prof, Tf0, controls, vialsize, fillvol

dom = Domain(config)
params, ncontrols = params_nondim_setup(cparams, controls)
um = LSS.make_u0_ndim(config)
ϕm = ϕ_T_from_u_view(um, dom)[1]
ϕm .+= .8*dom.rmax - 1e-6 + 1e-8

R = dom.rmax
L = dom.zmax
Ri = LSS.get_subf_r(ϕm, dom)

# LSS.eval_b(Tdm, Tdm, params)



p_num = solve_p(um, Tfm, Tdm, dom, params)
p_sol1 = gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=250)
p_anl = [LSS.calc_psub(ustrip(u"K", T0)) + (r < Ri ? 0 : p_sol1(r, z))
            for r in dom.rgrid, z in dom.zgrid]
pl1 = heat(p_num .- p_anl, dom)

function analytical_error_relmax(Nmax)
    p_sol1 = gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=Nmax)
    p_sol2 = gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=Nmax+1)
    p_anl1 = [LSS.calc_psub(ustrip(u"K", T0)) + (r < Ri ? 0 : p_sol1(r, z))
                for r in dom.rgrid, z in dom.zgrid]
    p_anl2 = [LSS.calc_psub(ustrip(u"K", T0)) + (r < Ri ? 0 : p_sol2(r, z))
                for r in dom.rgrid, z in dom.zgrid]

    return maximum((p_anl2 .- p_anl1) ./ p_anl1)
end

@time analytical_error_relmax(150)

function gen_topprof(er)
    um = LSS.make_u0_ndim(config)
    ϕm = ϕ_T_from_u_view(um, dom)[1]
    ϕm .+= .8dom.rmax + er - 0.99e-6#-2e-6

    R = dom.rmax
    L = dom.zmax
    Ri = LSS.get_subf_r(ϕm, dom)

    p_num = solve_p(um, Tfm, Tdm, dom, params)
    p_sol1 = gen_psol(Ri, R, L, Rp0, b, Δp; Nmax=150)
    p_anl = [LSS.calc_psub(ustrip(u"K", T0)) + (r < Ri ? 0 : p_sol1(r, z))
                for r in dom.rgrid, z in dom.zgrid]
    return [p_num[:,end], p_anl[:,end]]
end

perturbs = range(0, 2dom.dr, length=9) 
perturbs = perturbs .- step(perturbs)
@time res = [gen_topprof(ei) for ei in perturbs]
plot(res[1][1])
plot([ri[1][21] for ri in res])
plot!([ri[2][21] for ri in res])
plot(hcat(res...)')
res




pl1 = heat(p_num .- p_anl, dom)
pl2 = heat(p_num, dom)
pl3 = heat(p_anl, dom)

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

