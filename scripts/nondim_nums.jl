const LSS = LevelSetSublimation
include(scriptsdir("petr_exp_config.jl"))

b = eval_b(243.15u"K", LSS.calc_psub(243.15u"K"), config[:cparams])

αf = LSS.kf / LSS.ρ_ice / LSS.Cp_ice

αd = (LSS.k_sucrose*(1-LSS.ϵ_typical)) / (LSS.ρ_sucrose*(1-LSS.ϵ_typical)) / LSS.Cp_sucrose

R = 1u"cm"
L = 1u"cm"

t_thf = R^2/αf
t_thd = R^2/αd
t_m = b

t_sub = 5u"hr"

Fo_f = uconvert(NoUnits, t_sub/t_thf)
Fo_d = uconvert(NoUnits, t_sub/t_thd)
Fo_m = uconvert(NoUnits, t_sub/t_m)

@show Fo_f Fo_d Fo_m

ΔT = 40u"K"

m_p = 3u"g"
m_gl = 8u"g"

H_sens = (m_p*LSS.Cp_ice + m_gl*LSS.cp_gl) *ΔT
H_sub = LSS.ΔH * m_p

cakemass = (1-LSS.ϵ_typical)*LSS.ρ_sucrose*LSS.Cp_sucrose
ρvap = 1u"Torr"/LSS.R/233u"K"*LSS.Mw
vapmass =  ρvap* LSS.Cp_vap

Ste = uconvert(NoUnits, H_sens/H_sub)
thermalmass_ratio = uconvert(NoUnits, vapmass/cakemass)

flowthroughs = uconvert(NoUnits, LSS.ρ_ice/ρvap)

@show Ste thermalmass_ratio flowthroughs
