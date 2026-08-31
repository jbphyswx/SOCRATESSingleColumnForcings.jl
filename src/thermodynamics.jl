# ============================================================================================ #
# Thermodynamics
#
# The pipeline works directly on scalar fields `(T, p, q_tot)` — there is no thermodynamic state
# object. Derived quantities (condensate partition, density, virtual temperature, liquid-ice
# potential temperature) are recomputed on demand from those fields by a set of generic
# functions. The physics comes from a *backend* selected by dispatch on the `thermo_params`
# handle the caller threads through the pipeline:
#
#   * `DefaultThermodynamicsBackend` (defined here) — ideal-gas / Clausius–Clapeyron, usable with
#     no external thermodynamics dependency.
#   * a backend added by a package extension, which dispatches these functions on its own parameter set.
# ============================================================================================ #

"""
    AbstractThermodynamicsBackend

Supertype for the built-in, dependency-free thermodynamics backend. Extension backends dispatch
on their own parameter set and need not subtype this.
"""
abstract type AbstractThermodynamicsBackend end
Base.broadcastable(backend::AbstractThermodynamicsBackend) = tuple(backend)

"""
    DefaultThermodynamicsBackend()

Built-in dependency-free backend: ideal-gas density, Clausius–Clapeyron saturation (constant
`L_v0`/`L_s0`), and a bisection saturation adjustment. An extension backend can be used instead by
passing its parameter set as `thermodynamics_backend`.
"""
struct DefaultThermodynamicsBackend <: AbstractThermodynamicsBackend end


# --- generic methods: declared here, methods added per backend (default below; extension backends add theirs) --
"""Equilibrium condensate partition `(q_liq, q_ice)` from `(T, p, q_tot)`."""
function equilibrium_condensate end
"""Moist air density [kg/m³]."""
function air_density end
"""Virtual temperature [K]."""
function virtual_temperature end
"""Liquid-ice potential temperature [K]."""
function liquid_ice_pottemp end
"""Dry potential temperature [K]."""
function dry_pottemp end
"""Saturation adjustment from `(p, θ_liq_ice, q_tot)` → `(T, q_liq, q_ice)`."""
function saturation_adjust_pθq end
"""Saturation specific humidity over liquid."""
function q_vap_saturation_liq end
"""Saturation specific humidity over ice."""
function q_vap_saturation_ice end
"""Surface total specific humidity at saturation."""
function saturation_specific_humidity_from_pT end
"""Surface saturation total-water mixing ratio."""
function saturation_mixing_ratio_from_pT end
"""Saturation specific humidity from a saturation vapor pressure and a pressure."""
function q_vap_saturation_from_pressure end
"""Saturation specific humidity."""
function q_vap_saturation end
"""Saturation vapor pressure [Pa] for the phase whose latent heat at the triple point is given."""
function saturation_vapor_pressure end
"""Saturation vapor pressure over liquid water [Pa]."""
function saturation_vapor_pressure_liq end
"""Saturation vapor pressure over ice [Pa]."""
function saturation_vapor_pressure_ice end
"""Fraction of the condensate that is liquid."""
function liquid_fraction end
"""Dry-air gas constant [J/kg/K]."""
function R_d end
"""Water-vapor gas constant [J/kg/K]."""
function R_v end
"""Gravitational acceleration [m/s^2]."""
function grav end
"""Dry-air isobaric heat capacity [J/kg/K]."""
function cp_d end
"""Water-vapor isobaric heat capacity [J/kg/K]."""
function cp_v end
"""Liquid-water isobaric heat capacity [J/kg/K]."""
function cp_l end
"""Ice isobaric heat capacity [J/kg/K]."""
function cp_i end
"""Ratio of molar masses Mv/Md."""
function molmass_ratio end
"""Reference pressure for potential temperature [Pa]."""
function p_ref end
"""Freezing temperature [K]."""
function T_freeze end
"""All-ice threshold temperature [K]."""
function T_icenuc end
"""Latent heat of vaporization at the reference temperature [J/kg]."""
function L_v0 end
"""Latent heat of sublimation at the reference temperature [J/kg]."""
function L_s0 end
"""Saturation vapor pressure at the triple point [Pa]; the anchor shared by both phases."""
function e_ref end
"""Reference temperature the latent heats are anchored at [K]."""
function T_0 end
"""Latent heat at `T` [J/kg] from its reference value `LH_0` and the heat-capacity difference `Δcp`."""
function latent_heat_generic end
"""Latent heat of vaporization at `T` [J/kg]."""
function latent_heat_vapor end
"""Latent heat of sublimation at `T` [J/kg]."""
function latent_heat_sublim end

abstract type AbstractPhase end
struct Vapor <: AbstractPhase end
struct Liquid <: AbstractPhase end
struct Ice <: AbstractPhase end

# ============================================================================================ #
# Default backend physics (`DefaultThermodynamicsBackend`).
#
# Dependency-free ideal-gas / Clausius–Clapeyron approximations, used when no thermodynamics
# extension is loaded. A loaded extension backend defines these methods on its own parameter set.
# ============================================================================================ #

# Physical constants (SI), returned as `Float64`; functions convert to the working `FT`.
@inline e_ref(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(611.2)        # saturation vapor pressure at the triple point [Pa]
@inline R_d(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(287.058)       # dry-air gas constant [J/kg/K]
@inline R_v(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(461.5)               # water-vapor gas constant [J/kg/K]
@inline grav(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(9.80665)         # gravitational acceleration [m/s^2]
@inline cp_d(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(ATLAS_CP_D)          # dry-air heat capacity [J/kg/K]
@inline cp_v(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(1850.0)          # water-vapor heat capacity [J/kg/K]
@inline cp_l(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(4218.0)          # liquid-water heat capacity [J/kg/K]
@inline cp_i(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(2100.0)          # ice heat capacity [J/kg/K]
@inline molmass_ratio(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(0.622)  # M_v / M_d
@inline p_ref(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(1.0e5)          # reference pressure [Pa]
@inline T_0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(273.16)           # Reference temperatre [K]
@inline T_freeze(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(273.15)      # freezing temperature [K]
@inline T_icenuc(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(233.15)      # all-ice threshold [K] (−40 °C)
@inline L_v0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(ATLAS_LH_V0)     # latent heat of vaporization at T_freeze [J/kg]
@inline L_s0(::DefaultThermodynamicsBackend, ::Type{FT} = Float64) where {FT} = FT(ATLAS_LH_S0)     # latent heat of sublimation at T_freeze [J/kg]


@inline latent_heat_generic(backend::DefaultThermodynamicsBackend, T::FT, LH_0::FT, Δcp::FT) where {FT} =
    LH_0 + Δcp * (T - T_0(backend, FT))
@inline latent_heat_vapor(backend::DefaultThermodynamicsBackend, T::FT) where {FT} =
    latent_heat_generic(backend, T, L_v0(backend, FT), cp_v(backend, FT) - cp_l(backend, FT))
@inline latent_heat_sublim(backend::DefaultThermodynamicsBackend, T::FT) where {FT} =
    latent_heat_generic(backend, T, L_s0(backend, FT), cp_v(backend, FT) - cp_i(backend, FT))

"""
    saturation_vapor_pressure_liq(backend, T)

Saturation vapor pressure over liquid water [Pa], integrated Clausius–Clapeyron with constant
`L_v0(backend)`, anchored at `(T_freeze, 611.2 Pa)`.
"""
@inline function saturation_vapor_pressure_liq(backend::DefaultThermodynamicsBackend, T::FT;
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
) where {FT}
    return e_ref * exp(L_v0(backend, FT) / R_v(backend, FT) * (one(FT) / T_fr - one(FT) / T))
end

"""
    saturation_vapor_pressure_ice(backend, T)

Saturation vapor pressure over ice [Pa], integrated Clausius–Clapeyron with constant
`L_s0(backend)` (sublimation), anchored at `(T_freeze, 611.2 Pa)`.
"""
@inline function saturation_vapor_pressure_ice(backend::DefaultThermodynamicsBackend, T::FT;
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
) where {FT}
    return e_ref * exp(L_s0(backend, FT) / R_v(backend, FT) * (one(FT) / T_fr - one(FT) / T))
end

@inline function saturation_vapor_pressure(backend::DefaultThermodynamicsBackend, T::FT, phase::AbstractPhase = Liquid();
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
) where {FT}
    if phase === Liquid()
        return saturation_vapor_pressure_liq(backend, T; T_fr, e_ref)
    elseif phase === Ice()
        return saturation_vapor_pressure_ice(backend, T; T_fr, e_ref)
    else
        error("saturation_vapor_pressure: phase must be Liquid() or Ice()")
    end
end

"""
    q_vap_saturation_from_pressure(backend, e_sat, p)

Saturation specific humidity implied by a saturation vapor pressure `e_sat` at pressure `p`.
"""
@inline function q_vap_saturation_from_pressure(backend::DefaultThermodynamicsBackend, e_sat::FT, p::FT;
     ε::FT = molmass_ratio(backend, FT),
) where {FT}
    # At low pressure (high altitude) e_sat can approach/exceed p, driving the denominator ≤ 0. Physically
    # the air is then unsaturatable, so q* ≥ 1 (⟹ no condensate for any specific humidity ≤ 1).
    denom = p - (one(FT) - ε) * e_sat
    return denom > zero(FT) ? ε * e_sat / denom : one(FT)
end

"""
    q_vap_saturation_liq(backend, T, p)

Saturation specific humidity over liquid water.
"""
@inline function q_vap_saturation_liq(backend::DefaultThermodynamicsBackend, T::FT, p::FT;
     ε::FT = molmass_ratio(backend, FT),
     T_fr::FT = T_freeze(backend, FT),
     e_ref::FT = e_ref(backend, FT),
) where {FT}
    e_sat = saturation_vapor_pressure_liq(backend, T; T_fr, e_ref)
    return q_vap_saturation_from_pressure(backend, e_sat, p; ε)
end

"""
    q_vap_saturation_ice(backend, T, p)

Saturation specific humidity over ice.
"""
@inline function q_vap_saturation_ice(backend::DefaultThermodynamicsBackend, T::FT, p::FT;
     ε::FT = molmass_ratio(backend, FT),
     T_fr::FT = T_freeze(backend, FT),
     e_ref::FT = e_ref(backend, FT),
) where {FT}
    e_sat = saturation_vapor_pressure_ice(backend, T; T_fr, e_ref)
    return q_vap_saturation_from_pressure(backend, e_sat, p; ε)
end

"""
    q_vap_saturation(backend, T, p, phase)

Saturation specific humidity over a single `phase` ([`Liquid`](@ref) or [`Ice`](@ref)).
"""
@inline function q_vap_saturation(backend::DefaultThermodynamicsBackend, T::FT, p::FT, phase::AbstractPhase;
     ε::FT = molmass_ratio(backend, FT),
     T_fr::FT = T_freeze(backend, FT),
     e_ref::FT = e_ref(backend, FT),
) where {FT}
    e_sat = saturation_vapor_pressure(backend, T, phase; T_fr, e_ref)
    return q_vap_saturation_from_pressure(backend, e_sat, p; ε)
end



"""
    liquid_fraction(backend, T)

Fraction of condensate that is liquid: 1 at/above freezing, 0 at/below the all-ice threshold,
linear in temperature between.
"""
@inline function liquid_fraction(backend::DefaultThermodynamicsBackend, T::FT;
    T_fr::FT = T_freeze(backend, FT),
    T_icenuc::FT = T_icenuc(backend, FT),
) where {FT}
    return clamp((T - T_icenuc) / (T_fr - T_icenuc), zero(FT), one(FT))
end

"""
    q_vap_saturation(backend, T, p)

Phase-blended saturation specific humidity (liquid/ice weighted by [`liquid_fraction`](@ref)).
"""
@inline function q_vap_saturation(backend::DefaultThermodynamicsBackend, T::FT, p::FT;
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    T_icenuc::FT = T_icenuc(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    λ::FT = liquid_fraction(backend, T; T_fr, T_icenuc),
) where {FT}
    e_sat = λ * saturation_vapor_pressure_liq(backend, T; T_fr, e_ref) + (one(FT) - λ) * saturation_vapor_pressure_ice(backend, T; T_fr, e_ref)
    return q_vap_saturation_from_pressure(backend, e_sat, p; ε)
end

"""
    equilibrium_condensate(backend, T, p, q_tot) -> (; q_liq, q_ice)

Mixed-phase equilibrium condensate: total water above saturation condenses and is partitioned
into liquid and ice by [`liquid_fraction`](@ref).
"""
@inline function equilibrium_condensate(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT;
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    T_icenuc::FT = T_icenuc(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    λ::FT = liquid_fraction(backend, T; T_fr, T_icenuc),
) where {FT}
    q_c = max(zero(FT), q_tot - q_vap_saturation(backend, T, p; ε, T_fr, T_icenuc, e_ref, λ))
    return (; q_liq = λ * q_c, q_ice = (one(FT) - λ) * q_c)
end

@inline function virtual_temperature(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT;
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    (; q_liq, q_ice) = equilibrium_condensate(backend, T, p, q_tot; ε, T_fr, e_ref)
    return virtual_temperature(backend, T, q_tot, q_liq, q_ice; R_v, R_d)
end

@inline function virtual_temperature(backend::DefaultThermodynamicsBackend, T::FT, q_tot::FT, q_liq::FT, q_ice::FT;
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    q_vap = q_tot - q_liq - q_ice
    return T * (one(FT) + (R_v / R_d - one(FT)) * q_vap - q_liq - q_ice)
end

@inline function air_density(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT, q_liq::FT, q_ice::FT;
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    return p / (R_d * virtual_temperature(backend, T, q_tot, q_liq, q_ice; R_v, R_d))
end

@inline function air_density(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT;
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    R_v::FT = R_v(backend, FT),
    R_d::FT = R_d(backend, FT),
) where {FT}
    (; q_liq, q_ice) = equilibrium_condensate(backend, T, p, q_tot; ε, T_fr, e_ref)
    return air_density(backend, T, p, q_tot, q_liq, q_ice; R_v, R_d)
end

@inline function dry_pottemp(backend::DefaultThermodynamicsBackend, T::FT, p::FT;
    cp_d::FT = cp_d(backend, FT),
    p_ref::FT = p_ref(backend, FT),
) where {FT}
    κ = R_d(backend, FT) / cp_d
    return T * (p_ref / p)^κ
end

"""
    liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice)

Liquid-ice potential temperature: the composition-weighted moist potential temperature reduced
by the latent heat of the liquid and ice condensate (Tripoli–Cotton form).
"""
@inline function liquid_ice_pottemp(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT, q_liq::FT, q_ice::FT;
    cp_d::FT = cp_d(backend, FT),
    cp_v::FT = cp_v(backend, FT),
    cp_l::FT = cp_l(backend, FT),
    cp_i::FT = cp_i(backend, FT),
    p_ref::FT = p_ref(backend, FT),
    cp_m::FT = cp_d * (one(FT) - q_tot) + cp_v * (q_tot - q_liq - q_ice) + cp_l * q_liq + cp_i * q_ice,
    κ::FT = R_d(backend, FT) / cp_m,
) where {FT}
    return T * (p_ref / p)^κ * exp(-(L_v0(backend, FT) * q_liq + L_s0(backend, FT) * q_ice) / (cp_m * T))
end

function liquid_ice_pottemp(backend::DefaultThermodynamicsBackend, T::FT, p::FT, q_tot::FT;
    cp_d::FT = cp_d(backend, FT),
    cp_v::FT = cp_v(backend, FT),
    cp_l::FT = cp_l(backend, FT),
    cp_i::FT = cp_i(backend, FT),
    p_ref::FT = p_ref(backend, FT),
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    T_icenuc::FT = T_icenuc(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    λ::FT = liquid_fraction(backend, T; T_fr, T_icenuc),
) where {FT}
    (; q_liq, q_ice) = equilibrium_condensate(backend, T, p, q_tot; ε, T_fr, T_icenuc, e_ref, λ = λ)
    return liquid_ice_pottemp(backend, T, p, q_tot, q_liq, q_ice; cp_d, cp_v, cp_l, cp_i, p_ref)
end

"""
    saturation_adjust_pθq(backend, p, θ_liq_ice, q_tot) -> (T, q_liq, q_ice)

Saturation adjustment by bisection on temperature so that
`liquid_ice_pottemp(backend, T, p, q_tot) == θ_liq_ice`.
"""
function saturation_adjust_pθq(backend::DefaultThermodynamicsBackend, p::FT, θ_liq_ice::FT, q_tot::FT;
    maxiter::Int = 50, tol::FT = FT(1e-6),
    cp_d::FT = cp_d(backend, FT),
    cp_v::FT = cp_v(backend, FT),
    cp_l::FT = cp_l(backend, FT),
    cp_i::FT = cp_i(backend, FT),
    p_ref::FT = p_ref(backend, FT),
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    T_icenuc::FT = T_icenuc(backend, FT),
    e_ref::FT = e_ref(backend, FT),
    λ::FT = FT(NaN)
) where {FT}
    T_low, T_high = FT(150), FT(350)
    T = (T_low + T_high) / FT(2)
    for _ in 1:maxiter
        T = (T_low + T_high) / FT(2)
        λ = isnan(λ) ? liquid_fraction(backend, T; T_fr, T_icenuc) : λ
        θ_mid = liquid_ice_pottemp(backend, T, p, q_tot; cp_d, cp_v, cp_l, cp_i, p_ref, ε, T_fr, T_icenuc, e_ref, λ = λ)
        (abs(θ_mid - θ_liq_ice) < tol) && break
        θ_mid < θ_liq_ice ? (T_low = T) : (T_high = T)
    end
    (; q_liq, q_ice) = equilibrium_condensate(backend, T, p, q_tot; ε, T_fr, T_icenuc, e_ref, λ = λ)
    return (; T, q_liq, q_ice)
end

"""
    saturation_specific_humidity_from_pT(backend, pg, Tg)

Surface total specific humidity at saturation over liquid, given surface pressure `pg` and
temperature `Tg`.
"""
@inline function saturation_specific_humidity_from_pT(backend::DefaultThermodynamicsBackend, p::FT, T::FT, phase::AbstractPhase = Liquid();
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
) where {FT}
    q_vap_saturation(backend, T, p, phase; ε, T_fr, e_ref)
end

"""
    saturation_mixing_ratio_from_pT(backend, p, T)

Surface saturation total-water content at saturation over liquid:
`w_s = ε · e_sat / (pg − e_sat)`, with `ε = molmass_ratio = M_v/M_d ≈ 0.622`.
"""
@inline function saturation_mixing_ratio_from_pT(backend::DefaultThermodynamicsBackend, p::FT, T::FT, phase::AbstractPhase = Liquid();
    ε::FT = molmass_ratio(backend, FT),
    T_fr::FT = T_freeze(backend, FT),
    e_ref::FT = e_ref(backend, FT),
) where {FT}
    pv = saturation_vapor_pressure(backend, T, phase; T_fr, e_ref)
    return ε * pv / (p - pv)
end
