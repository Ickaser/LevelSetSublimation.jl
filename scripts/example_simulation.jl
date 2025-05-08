

# -------------------------
# Some base physical properties, which only change if using
# vials which are not borosilicate or if using solvent other than water
base_props = PhysicalProperties()

# ----------------------------
# Simulation geometry inputs: these names must be exact
vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (41, 31)

dom = Domain(@dict vialsize fillvol simgridsize)

# -------------------------
# Case-specific inputs, documented in TimeConstantProperties
# These names can be changed freely, as long as they are passed
# to the TimeConstantProperties constructor in the right order

Rp0 = 1.4u"cm^2*hr*Torr/g"

# Varying pore structure in vertical direction
A1 = 16u"cm*hr*Torr/g"
A2 = 0.5u"1/cm"
Tguess = 260u"K"
hd = (dom.zmax .- dom.zgrid)*u"m"
dRpdh = @. A1/(1 + A2*hd)^2 # This is a constant for A2=0
# If l is a constant everywhere, it can be passed as a single value rather than array.
# Here we vary the pore structure, so it as an array of values
l = [sqrt(base_props.R*Tguess/base_props.Mw) / dRp for r in dom.rgrid, dRp in dRpdh] # Varying in height

c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
ϵ = (ρ_solution-c_solid)/ρ_solution 
κ = 0.0u"m^2" # no viscous-regime flow allowed
kd = LevelSetSublimation.k_sucrose*(1-ϵ)
Kvwf = 20.0u"W/m^2/K"
m_v = get_vial_mass(vialsize)
A_v = get_vial_radii(vialsize)[2]^2 *π
B_d = 0u"Ω/m^2"
B_f = 5e7u"Ω/m^2"
B_vw = 2e6u"Ω/m^2"
# Assemble
tcp = TimeConstantProperties(
    ϵ, l, κ, Rp0, kd, Kvwf, 
    m_v, A_v, B_d, B_f, B_vw)

# ------------------------------
# Process parameters which may change over time.
# Documented in TimeVaryingProperties; names can be chosen, but order must be matched
# RampedVariable is a convenience type from LyoPronto.jl, documented there
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 1u"K/minute")
pch = RampedVariable(100u"mTorr")
f_RF = RampedVariable(12u"GHz")
P_per_vial = RampedVariable(0.5u"W")
# Kv
Kc = 3.58e-4u"cal/s/K/cm^2" 
Kd = 11.6e-4u"cal/s/K/cm^2/Torr"
Kp = 0.46u"1/Torr"
Kshf = RpFormFit(Kc, Kd, Kp)
# Assemble
tvp = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

# All parameters put together in one tuple
paramsd = (base_props, tcp, tvp)

# ---------------------------------
# Pack everything into a dictionary for simulating or saving
# The dictionary keys need to be exact, which is
# This is one way of constructing the dictionary, if your variable 
# names match the required dictionary keys
config = @dict paramsd vialsize fillvol simgridsize
config[:time_integ] = :dae_then_exp

alt_simgridsize = (41, 31)
# Or we can do the following, if variables have other names
config = Dict{Symbol, Any}(
    :paramsd => paramsd,
    :vialsize => vialsize,
    :fillvol => fillvol,
    :simgridsize => alt_simgridsize,
    :time_integ => :dae_then_exp
)

# -----------------------------
# Now, to run the simulation:

# To just run it once and see progress in stdout...
res = sim_from_dict(config, verbose=true)

# or, check first to see if a simulation with this config has been run,
# loading those results if it has already been run.
# if not already run, simulation runs quietly
# See DrWatson documentation for more `produce_or_load`
res, fname = produce_or_load(sim_from_dict, config; 
    filename=hash, prefix=datadir("sims", "ex"))

sim = res["sim"]
# -------------------------------
# Generate some plots

# First, plot temperature history for several space points
locs = [(0.0, 0.0), (0.2, 0.4), (0.7, 0.7)]
vtmarks = [:diamond, :utriangle, :heptagon, :dtriangle]
blankplothrC()
vt_plot!(sim, locs, markers=permutedims(vtmarks))

# Then, show the temperature field at several time points
summaryT(sim)
# Add some markers to that plot showing 
# where the temperature histories came from
placethermocouples!(dom, locs, c=palette(:Oranges_4)[4:-1:2],
    markers=vtmarks, markersize=8)
