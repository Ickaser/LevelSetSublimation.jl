const LSS = LevelSetSublimation

using SpecialFunctions
using Roots
using Unitful
using Plots
using FastGaussQuadrature

R = .01 # 1cm
Ri = .005 #0.5cm
L = .02 #2cm
Rp0 = ustrip(u"m/s", 1.4u"Torr*cm^2/g*hr")
A1 = 16u"cm*hr*Torr/g"
Mw = 18u"g/mol"
Rg = 8.3145u"J/mol/K"
T0 = 250u"K"
l_bulk = sqrt(Rg*T0/Mw) / A1
b = ustrip(u"s", l_bulk*sqrt(Mw/Rg/T0)) #Mw/R/T * (l*NaNMath.sqrt(R*T/Mw)
Δp = ustrip(u"mTorr", LSS.calc_psub(T0) - 100u"mTorr") # Subtract off chamber pressure

Nmax = 50
charfun_r(l) = besselj0(l*Ri)*bessely1(l*R) - bessely0(l*Ri)*besselj1(l*R)
λj = [find_zero(charfun_r, (2j-1)*π/R) for j in 1:Nmax]

charfun_z(m) = cos(m*L) - Rp0*b*m*sin(m*L)
@show π/L
μk = permutedims([find_zero(charfun_z, (k-0.5)*π/L) for k in 1:Nmax])

xx = range(1, 500)
plot(xx, charfun_z.(xx))
	
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


x, w = gausslegendre(50)
r_nodes = @. (x +1)/2*(R-Ri) + Ri
r_wts = @. w/2*(R-Ri)
z_nodes = @. (x +1)/2*L
z_wts = @. w/2*L

innerprod_r(f, g) = sum(@. r_wts * r_nodes * f(r_nodes) * g(r_nodes))
innerprod_z(f, g) = sum(@. z_wts * f(z_nodes) * g(z_nodes))
u_normsq = [innerprod_r(uj, uj) for uj in u]
v_normsq = [innerprod_z(vk, vk) for vk in v]

j = 1:Nmax
k = permutedims(1:Nmax)

# # Using p* = p-pch
# Cj = [Ri*make_duj(jj)(Ri)*Δp for jj in j]
# # Cjk = @. Cj*sin(μk*L)/μk
# Cjk = [innerprod_z(x->Ci, vk) for Ci in Cj, vk in v]
# Chjk = @. Cjk / (λj^2 + μk^2)
# coeffs = Chjk .* u_normsq .* permutedims(v_normsq)

# Using p' = p-psub, r-transform then z-transform
# Ck = permutedims([vk(L) for vk in v]) .* (-Δp/Rp0/b)
# Cjk = @. Ck  / (λj^2 + μk^2)
# coeffs = Cjk ./ u_normsq ./ permutedims(v_normsq)

# Using p' = p-psub, z-transform then r-transform
Ck = [vk(L) for vk in v] .* (-Δp/Rp0/b)
Cjk = [innerprod_r(x->Ci, uj) for uj in u, Ci in Ck] 
coeffs = Cjk ./ u_normsq ./ permutedims(v_normsq)


[uj(R) for uj in u]
[uj(Ri) for uj in u]
[vk(0) for vk in v]
[vk(L) for vk in v]

function p_sol(r, z)
    sum(coeffs .* [uj(r) for uj in u] .* permutedims([vk(z) for vk in v]))
end
@time p_sol(R*0.75, L/2)

rr = range(0.005, .01, length=100)
zz = range(0, L, length=110)
@time ps = [p_sol(r, z) for z in zz, r in rr]
heatmap(rr, zz, ps)
plot(rr, [uj.(rr) for uj in u])
plot([vk.(zz) for vk in v])
ps[:,1]
