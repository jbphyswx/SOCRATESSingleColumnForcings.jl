module SOCRATESSingleColumnForcingsThermodynamicsExt

# Thermodynamics.jl extension backend: Thermodynamics.jl
# Each method builds a Thermodynamics state (`PhaseEquil_pTq` / `PhaseEquil_pθq`)
# or a `PhasePartition` and queries it.  The default fallback backend is in core.

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD

const TD_version = pkgversion(TD)

# --- physical constants ---
@inline SSCF.R_d(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params))
@inline SSCF.R_v(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_v(thermo_params))
@inline SSCF.grav(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.grav(thermo_params))
# ε = M_v/M_d ≈ 0.622 = R_d/R_v. Thermodynamics' own `molmass_ratio` is the INVERSE (M_d/M_v ≈ 1.608),
# so compute ε directly from R_d/R_v — matches the backend contract and is version-independent.
@inline SSCF.molmass_ratio(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params) / TD.Parameters.R_v(thermo_params))

# --- equilibrium condensate partition from (T, p, q_tot) ---
function SSCF.equilibrium_condensate(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    q = TD.PhasePartition(thermo_params, TD.PhaseEquil_pTq(thermo_params, p, T, q_tot))
    return (; q_liq = q.liq, q_ice = q.ice)
end

# --- air density ---
# condensate-aware (given partition): ρ = p / (R_m(q)·T)
SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice) =
    TD.air_density(thermo_params, T, p, TD.PhasePartition(q_tot, q_liq, q_ice))
# from total water (equilibrium partition resolved internally)
SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot) =
    TD.air_density(thermo_params, TD.PhaseEquil_pTq(thermo_params, p, T, q_tot))

# --- virtual temperature ---
function SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.virtual_temperature(thermo_params, T, q_tot, q_liq, q_ice)
end

if TD_version ≥ v"0.11.7" # https://github.com/CliMA/Thermodynamics.jl/pull/191/changes
    # (T, q_tot, q_liq, q_ice) -> T_v; 0.11 needs no ρ (T_v = R_m(q)/R_d · T)
    SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, q_tot, q_liq, q_ice) =
        TD.virtual_temperature(thermo_params, T, TD.PhasePartition(q_tot, q_liq, q_ice))
else
    # (T, q_tot, q_liq, q_ice) -> T_v; < 0.11.7 takes an (unused) positional ρ (T_v = R_m(q)/R_d · T)
    SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, q_tot, q_liq, q_ice) =
        TD.virtual_temperature(thermo_params, T, zero(T), TD.PhasePartition(q_tot, q_liq, q_ice))
end

# --- liquid-ice potential temperature (from pressure + resolved partition) ---
SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice) =
    TD.liquid_ice_pottemp_given_pressure(thermo_params, T, p, TD.PhasePartition(q_tot, q_liq, q_ice))
function SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.liquid_ice_pottemp(thermo_params, T, p, q_tot, q_liq, q_ice)
end

# --- dry potential temperature ---
# `dry_pottemp_given_pressure` was renamed to `potential_temperature_given_pressure` in the Thermo 0.15.3
# restructure (PR #293); both are the dry (default-q) potential temperature from pressure.
if TD_version < v"0.15.3"
    SSCF.dry_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) = TD.dry_pottemp_given_pressure(thermo_params, T, p)
else
    SSCF.dry_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) = TD.potential_temperature_given_pressure(thermo_params, T, p)
end

# --- saturation adjustment (p, θ_liq_ice, q_tot) -> (; T, q_liq, q_ice) ---
# `PhaseEquil_pθq` runs the θ_liq_ice saturation adjustment internally. `tol === nothing` lets
# Thermodynamics use its own default relative tolerance.
function SSCF.saturation_adjust_pθq(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, p, θ_liq_ice, q_tot; maxiter::Int = 50, tol = nothing) where {FT}
    _tol = tol === nothing ? nothing : FT(tol)
    ts = TD.PhaseEquil_pθq(thermo_params, FT(p), FT(θ_liq_ice), FT(q_tot), maxiter, _tol)
    q = TD.PhasePartition(thermo_params, ts)
    return (; T = TD.air_temperature(thermo_params, ts), q_liq = q.liq, q_ice = q.ice)
end

# --- saturation specific humidity over liquid ---
# `q_vap_saturation_generic` was renamed to `q_vap_saturation(…, phase)` in the Thermo 0.15.3 restructure (#293).
if TD_version < v"0.15.3"
    _q_vap_saturation_liquid(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ρ) = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
else
    _q_vap_saturation_liquid(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ρ) = TD.q_vap_saturation(thermo_params, T, ρ, TD.Liquid())
end

SSCF.q_vap_saturation_liquid(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    _q_vap_saturation_liquid(thermo_params, T, TD.air_density(thermo_params, T, p))

function SSCF.calc_qg_from_pgTg(thermo_params::TD.Parameters.ThermodynamicsParameters, pg, Tg)
    ρg = TD.air_density(thermo_params, Tg, pg)  # dry surface density
    return _q_vap_saturation_liquid(thermo_params, Tg, ρg)
end

# --- surface saturation total-water mixing ratio: w_s = ε·e_s/(pg − e_s), ε = M_v/M_d ---
function SSCF.saturation_q_tot_from_pgTg(thermo_params::TD.Parameters.ThermodynamicsParameters, pg, Tg)
    pvg = TD.saturation_vapor_pressure(thermo_params, Tg, TD.Liquid())
    ε = SSCF.molmass_ratio(thermo_params)
    return ε * pvg / (pg - pvg)
end

end # module
