# dom = Domain(51, 51, 1.0, 1.0)
cparams = p = make_default_params()
ϕ0type = :circ
Tf0 = 233.15u"K"
Q_gl_RF = 0.002u"W" # = volumetric * relevant vial volume
t_samp = (0:0.1:1) .* u"hr"
# T_sh = [233.15
Tsh = 263.15u"K"
Q_ic = 0.3u"W/cm^3"
p_ch = 100u"mTorr"

controls = Dict{Symbol, Any}()
@pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch


vialsize = "10R"
fillvol = 2u"mL"

config = Dict{Symbol, Any}()
@pack! config = cparams, ϕ0type, Tf0, controls, vialsize, fillvol


# Set up stuff to make debugging easier
params, meas_keys, ncontrols = params_nondim_setup(cparams, controls)

r_vial = get_vial_rad(vialsize)
z_fill = fillvol / π / r_vial^2

rmax = ustrip(u"m", r_vial)
zmax = ustrip(u"m", z_fill)

dom = Domain(51, 51, rmax, zmax)

ϕ0 = make_ϕ0(ϕ0type, dom)   
u0 = zeros(dom.ntot+2)
u0[1:dom.ntot] = reshape(ϕ0, :)
u0[dom.ntot+1] = ustrip(u"K", Tf0)
u0[dom.ntot+2] = ustrip(u"K", Tf0)
reinitialize_ϕ_HCR!(ϕ0, dom)
T0 = solve_T(u0, dom, params)
p0 = solve_p(u0, T0, dom, params)
