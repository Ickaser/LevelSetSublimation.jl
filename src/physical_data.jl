
# A listing of physical properties for use in simulation, with references to sources

# -----------------------------------
# Sublimation pressure, to more significant figures than necessary but hey
# Left without units becuase this is used inside functions


# Triple point, without units: from IAPWS 1995, as rereleased in 2018
const psub_Tt = 273.16 # Triple point T, kelvin
const psub_pt = 611.654_771 # Triple point p, Pa


const psub_ei = [
 20.996_966_510_789_7,
  3.724_374_782_713_62,
-13.920_548_321_552_4,
 29.698_876_501_356_6,
-40.197_239_263_594_4,
 29.788_048_105_021_5,
 -9.130_509_635_477_21,
]
"""
    calc_psub_highorder(T)

Compute pressure (in Pascals) of sublimation at temperature `T` in Kelvin

From Feistel and Wagner, 2006
"""
function calc_psub_highorder(T)
    # if clamp(T, 20.0, 273.0) != T
    #     @warn "Invalid temperature for sublimation pressure correlation. Clamped to 20 and 273 K" T
    #     T = clamp(T, 20.0, 273.0)
    # end
    η = sum(psub_ei .* (T/psub_Tt) .^ (0:6))
    lnπ = 1.5log(T/psub_Tt) + (1 - psub_Tt/T) * η
    return exp(lnπ) * psub_pt
end

const R = Unitful.R
import LyoPronto: ρ_ice, cp_ice, Mw, μ_vap, k_ice, ΔH, calc_psub
const μ = LyoPronto.μ_vap
const kf = LyoPronto.k_ice
import LyoPronto: ρ_gl, k_gl, cp_gl, εpp_gl, e_0
const ε0 = e_0
import LyoPronto: ρ_sucrose, k_sucrose

const ρ_wat = 997.0u"kg/m^3"

# Water vapor heat capacity near triple point
const Cp_vap = 1.86u"kJ/kg/K"

# From NIST, multiple sources
const Cp_sucrose = 422.50u"J/mol/K" / 342.3u"g/mol"

# From a solar paper and nanomaterials paper.
# Not clear what temperature (or even state) this is for
const k_mannitol1 = 0.60u"W/m/K"
const k_mannitol2 = 0.74u"W/m/K"

# Typical value of porosity, heuristically
# const ϵ_typical = 0.92

# Guessed values of porous structure parameters
const l_guess = 1e-5 * u"m" # 10 micron pore size
const κ_guess = 1e-10 * u"m^2" # 100 Darcys

# Typical value at around pch = 20 Pa
const Kshf_guess = 15.0u"W/m^2/K"
const pch_guess = 20.0u"Pa"

