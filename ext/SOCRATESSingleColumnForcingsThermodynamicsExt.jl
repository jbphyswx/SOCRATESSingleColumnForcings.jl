module SOCRATESSingleColumnForcingsThermodynamicsExt

# Thermodynamics.jl backend. The SSCF thermodynamics methods defined here dispatch on
# `ThermodynamicsParameters`, so the parameter set itself serves as the `thermodynamics_backend`:
# plain scalar-field calls into Thermodynamics.jl, no state object. The built-in default backend
# `DefaultThermodynamicsBackend` lives in the core.

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD

# --- physical constants ---
@inline SSCF.R_d(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params))
@inline SSCF.R_v(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_v(thermo_params))
@inline SSCF.grav(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.grav(thermo_params))
# SSCF's default uses ε = M_v/M_d ≈ 0.622 (meteorological), e.g. q_sat = ε·e/(p−(1−ε)e). Thermodynamics'
# own `molmass_ratio` is the INVERSE (M_d/M_v = R_v/R_d ≈ 1.608, aliased to `Rv_over_Rd`). Compute ε
# directly from R_d/R_v so the accessor matches the backend contract regardless of API version.
@inline SSCF.molmass_ratio(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params) / TD.Parameters.R_v(thermo_params))

# --- equilibrium condensate partition from (T, p, q_tot) ---
function SSCF.equilibrium_condensate(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    ρ = TD.air_density(thermo_params, T, p, q_tot)
    (q_liq, q_ice) = TD.condensate_partition(thermo_params, T, ρ, q_tot)  # (q_liq, q_ice)
    return (; q_liq, q_ice)
end

# ------------------------------------------------------------ #

function SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice)
    return TD.air_density(thermo_params, T, p, q_tot, q_liq, q_ice)
end

# (backend, T, p, q_tot) -> ρ — moist density from total water (the form the pipeline calls).
SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot) =
    TD.air_density(thermo_params, T, p, q_tot)

function SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.virtual_temperature(thermo_params, T, q_tot, q_liq, q_ice)
end

# ------------------------------------------------------------ #

# (backend, T, q_tot, q_liq, q_ice) -> T_v — matches the default backend's signature (no `p`; the
# partition suffices). `lev_to_z` calls this 4-argument (after-backend) form.
function SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, q_tot, q_liq, q_ice)
    return TD.virtual_temperature(thermo_params, T, q_tot, q_liq, q_ice)
end

# ------------------------------------------------------------ #

# liquid-ice potential temperature, computed from pressure with the resolved condensate partition
# (`liquid_ice_pottemp_given_pressure` takes `p`, not density).
function SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice)
    return TD.liquid_ice_pottemp_given_pressure(thermo_params, T, p, q_tot, q_liq, q_ice)
end
function SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.liquid_ice_pottemp(thermo_params, T, p, q_tot, q_liq, q_ice)
end

# ------------------------------------------------------------ #

SSCF.dry_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) = TD.potential_temperature_given_pressure(thermo_params, T, p)

# ------------------------------------------------------------ #


function SSCF.saturation_adjust_pθq(thermo_params::TD.Parameters.ThermodynamicsParameters, p, θ_liq_ice, q_tot; maxiter::Int = 50, tol = nothing)
    FT = eltype(thermo_params)
    _tol = tol === nothing ? FT(1.0e-6) : FT(tol)
    sat = TD.saturation_adjustment(TD.RS.NewtonsMethod, thermo_params, TD.pθ_li(), FT(p), FT(θ_liq_ice), FT(q_tot), maxiter, _tol)
    return (; T = sat.T, q_liq = sat.q_liq, q_ice = sat.q_ice)
end

SSCF.q_vap_saturation_liquid(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    TD.q_vap_saturation(thermo_params, T, TD.air_density(thermo_params, T, p), TD.Liquid())

function SSCF.calc_qg_from_pgTg(thermo_params::TD.Parameters.ThermodynamicsParameters, pg, Tg)
    ρg = TD.air_density(thermo_params, Tg, pg)
    return TD.q_vap_saturation(thermo_params, Tg, ρg, TD.Liquid())  # surface specific humidity over liquid
end

function SSCF.saturation_q_tot_from_pgTg(thermo_params::TD.Parameters.ThermodynamicsParameters, pg, Tg)
    pvg = TD.saturation_vapor_pressure(thermo_params, Tg, TD.Liquid())
    ε = SSCF.molmass_ratio(thermo_params)  # M_v/M_d ≈ 0.622 (see the accessor above; NOT TD's inverse)
    return ε * pvg / (pg - pvg)  # saturation total-water mixing ratio at the surface: w_s = ε·e_s/(p−e_s)
end

end # module
