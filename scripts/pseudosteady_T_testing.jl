# using BoundaryValueDiffEq
using SpecialFunctions
using FastGaussQuadrature
const LSS = LevelSetSublimation

# ------------ Simulation conditions


Rp0 = ustrip(u"m/s", 1.4u"Torr*cm^2/g*hr")
A1 = 16u"cm*hr*Torr/g"
Mw = 18u"g/mol"
Rg = 8.3145u"J/mol/K"
T0 = 250.15u"K"
l_bulk = sqrt(Rg*T0/Mw) / A1
b = ustrip(u"s", l_bulk*sqrt(Mw/Rg/T0)) #Mw/R/T * (l*NaNMath.sqrt(R*T/Mw)
Δp = ustrip(u"Pa", LSS.calc_psub(T0) - 100u"mTorr")
##
cparams = make_default_params()
cparams[:l] = upreferred(l_bulk)
cparams[:κ] = 0.0u"m^2"
cparams[:Rp0] = Rp0*u"m/s"
cparams[:Kv] *= 0
cparams[:Kw] = 10u"W/m^2/K"
##
simgridsize = (101, 101)
Tfm = fill(ustrip(u"K", T0), simgridsize[1])
Tdm = fill(ustrip(u"K", T0), simgridsize)
Tf0 = T0
Tw0 = T0 + 20u"K"
##
controls = Dict{Symbol, Any}()
Q_gl_RF = RampedVariable(0.0u"W")
Tsh = RampedVariable(0.0u"K")
Q_ic = RampedVariable(0.0u"W/cm^3")
p_ch = RampedVariable(100u"mTorr")
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch
### 
init_prof = :rad
vialsize = "6R"
fillvol = 5u"mL"
config = Dict{Symbol, Any}()
@pack! config = simgridsize, cparams, init_prof, Tf0, Tw0, controls, vialsize, fillvol
###
dom = Domain(config)
params, ncontrols = params_nondim_setup(cparams, controls)
um = LSS.make_u0_ndim(config)
ϕm = ϕ_T_from_u_view(um, dom)[1]
ϕm .+= .4*dom.rmax - 1e-6 +1e-8
###
R = dom.rmax
L = dom.zmax
Ri = LSS.get_subf_r(ϕm, dom)

# Tfs = @time LSS.pseudosteady_Tf(um, dom, params, fill(250.0, dom.nr))
# LSS.solve_T(um, Tfs, dom, params)

# --------------------------- Frozen T ODE

function ice_T_ode!(du, u, pars, r)
    T, dT = u
    ΔH, R0, p_ch, kf = pars
    # du[2] = 
    dTr1 = r > 0 ? dT/r : 0.0

    du[1] = dT
    du[2] = -dTr1 + ΔH/R0*(calc_psub(T) - p_ch)/kf/L
    # if any(isnan.(du))
    #     @info "problem" du u
    # end
end
# function ice_T_bc!(resid, u, pars, r)
#     Qpps = pars[4]
#     resid[1] = u[1][2] # dTdr = 0 at r=0
#     resid[2] = u[end][2] - Qpps # dTdr = - Q at r=0
# end

rspan = (0.0, Ri)

pars = (params[:ΔH], Rp0, ustrip(u"Pa", p_ch(0)), params[:kf])

u0 = [250, 0]

ode1 = ODEProblem(ice_T_ode!, u0, rspan)
solve(ode1, Tsit5(), p=pars, u0=[238.5, 0])

manual_shoot(T0, Qpps) = solve(ode1, Tsit5(), p=pars, u0=[T0, 0])[2,end] - Qpps
# bvp1 = TwoPointBVProblem(ice_T_ode!, ice_T_bc!, u0, rspan)
# sol1 = solve(bvp1, Shooting(Tsit5()), p=pars)

using Roots
sol1 = find_zero(x->manual_shoot(x, 30), 250) 



# -------------- Pressure gradient term

function gen_psol(Ri, R, L, Rp0, b, Δp; Nr=150, Nz=1000)
    charfun_r(l) = besselj0(l*Ri)*bessely1(l*R) - bessely0(l*Ri)*besselj1(l*R)
    λj = [find_zero(charfun_r, (j-0.5)*π/(R-Ri)) for j in 1:Nr]
    charfun_z(m) = cos(m*L) - Rp0*b*m*sin(m*L)
    # μk = permutedims([find_zero(charfun_z, (k-0.5)*π/L) for k in 1:Nmax])
    μk = zeros(1, Nz)
    μk[1] = find_zero(charfun_z, 0.5*π/L)
    for i in 2:Nz
        prev_root = μk[i-1]
        root1 = find_zero(charfun_z, prev_root + π/L, maxiters=100)
        if isapprox(root1, prev_root+π/L, rtol=.3)
            μk[i] = root1
        else
            @warn "problems" prev_root root1 prev_root+π/L
        end
    end
    # Evaluate eigenfunctions
    function make_uj(j)
        uj(r) = -bessely0(λj[j]*Ri)*besselj0(λj[j]*r) + besselj0(λj[j]*Ri)*bessely0(λj[j]*r)
    end
    function make_duj(j)
        duj(r) = λj[j]* bessely0(λj[j]*Ri)*besselj1(λj[j]*r) - λj[j]*besselj0(λj[j]*Ri)*bessely1(λj[j]*r)
    end
    u = [make_uj(j) for j in 1:Nr]
    du = [make_duj(j) for j in 1:Nr]
    function make_vk(k)
        vk(z) = cos(μk[k]*z)
    end
    function make_dvk(k)
        dvk(z) = -μk[k]*sin(μk[k]*z)
    end
    v = [make_vk(k) for k in 1:Nz]
    dv = [make_dvk(k) for k in 1:Nz]
    # Set up inner products and evaluate
    x, w = gausslegendre(Nr)
    r_nodes = @. (x +1)/2*(R-Ri) + Ri
    r_wts = @. w/2*(R-Ri)
    x, w = gausslegendre(Nz)
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

    function p_sol(r, z; coef = coeffs)
        sum(coef .* [uj(r) for uj in u] .* permutedims([vk(z) for vk in v]))
    end
    function dpdr_anl(r, z; coef = coeffs)
        sum(coef .* [duj(r) for duj in du] .* permutedims([vk(z) for vk in v]))
    end
    function dpdz_anl(r, z; coef = coeffs)
        sum(coef .* [uj(r) for uj in u] .* permutedims([dvk(z) for dvk in dv]))
    end
    return p_sol, dpdr_anl, dpdz_anl, coeffs
end

psol1, dpdr_1, dpdz_1, coeffs_p_anl = gen_psol(Ri, R, L, Rp0, b, Δp, Nr=300, Nz = 1500)
rescaled_coeffs(Δp_new, Δp_old) = coeffs_p_anl .* (Δp_new / Δp_old)

# Verify this method doesn't screw things up
psol1(dom.rmax, dom.zmax) ≈ psol1(dom.rmax, dom.zmax; coef = rescaled_coeffs(Δp, Δp))
dpdr_1(Ri, L/2) ≈ dpdr_1(Ri, L/2; coef = rescaled_coeffs(Δp, Δp))


# Sanity checks: mass conservation

N_pts = 100
x, w = gausslegendre(N_pts)
r_nodes = @. (x +1)/2*(R-Ri) + Ri
r_wts = @. w/2*(R-Ri)
z_nodes = @. (x +1)/2*L
z_wts = @. w/2*L
r_fluxes = [b*dpdr_1(Ri, zi) for zi in z_nodes]
tot_flow_anl = sum(r_fluxes .*z_wts)*Ri*2π
z_fluxes = [-(psol1(ri, dom.zmax) +Δp)/Rp0 for ri in r_nodes]
tot_flow_anl2 = sum(z_fluxes .*r_nodes .*r_wts ) *2π
(tot_flow_anl - tot_flow_anl2) ./ tot_flow_anl

Ωr = dom.rgrid .> Ri
tot_flow_num = sum(dom.rgrid[Ωr].*-(p_num[Ωr, end] .- ustrip(p_ch(0)))/Rp0/b)*dom.dr*b * 2π
tot_flow_num2 = sum((p_num[22, :] .- p_num[21, :])/dom.dr)*dom.dz*Ri*2π*b
# tot_flow_num2 = sum(LSS.compute_pderiv(um, )/dom.dr)*dom.dz*Ri*2π*b
(tot_flow_num - tot_flow_num2) ./ tot_flow_num

# Function for evaluating mass flux

function interface_mflux(Tsub)
    psub = LSS.calc_psub(Tsub)
    Δp_new = psub - ustrip(u"Pa", p_ch(0))
    r_fluxes = [b*dpdr_1(Ri, zi, coef = rescaled_coeffs(Δp_new, Δp)) for zi in z_nodes]
    # @info "mflux" extrema(r_fluxes)./b
    tot_flow_anl = sum(r_fluxes .*z_wts)*Ri*2π
    flux = tot_flow_anl / L / (Ri*2π)
end
    

# --- Analytical T 

Bi = uconvert(NoUnits, cparams[:Kw]*R*u"m"/cparams[:k])
function interface_Tflux(Tsub)
    C1 = Bi*(ustrip(u"K", Tw0)-Tsub)/(1+Bi*log(R/Ri))
    return params[:k]*C1/Ri
end
Tsub = 238.5
function analytical_T(r, Tsub)
    C1 = Bi*(ustrip(u"K", Tw0)-Tsub)/(1+Bi*log(R/Ri))
    return C1*log(r/Ri) + Tsub
end

interface_Tflux(238.5)
interface_mflux(238.5)*params[:ΔH]

nlobj(238.5)

analytical_T.(dom.rgrid[Ωr], Tsub)

# ----- Assemble into a nonlinear solve function for center T

function nlobj(T_c)
    Tf_sol = solve(ode1, Tsit5(), p=pars, u0=[T_c, 0])
    heat_f = - params[:kf] * Tf_sol[2,end] 
    heat_d = interface_Tflux(Tf_sol[1,end])
    mass_d = interface_mflux(Tf_sol[1,end])
    heat_l = params[:ΔH] * mass_d
    # @info "analytical" Tf_sol[1,end] Tf_sol[2,end] mass_d/b heat_d/params[:k] heat_f heat_d heat_l
    heat_f + heat_d + heat_l
end

nlobj(230)
@time Tc_sol = find_zero(nlobj, 230)

Tf_anl = solve(ode1, u0=[Tc_sol, 0], Tsit5(), p=pars)

Tfs = @time LSS.pseudosteady_Tf(um, dom, params, fill(250.0, dom.nr))
LSS.solve_T(um, Tfs, dom, params)[:,1]

Tfs[.~Ωr] .- Tf_anl.(dom.rgrid[.~Ωr], idxs=1)
Tf_anl.(dom.rgrid[.~Ωr], idxs=1)+0.5dom.dr
Tfs[.~Ωr]








