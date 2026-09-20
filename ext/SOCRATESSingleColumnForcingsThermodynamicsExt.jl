module SOCRATESSingleColumnForcingsThermodynamicsExt

# Thermodynamics.jl backend. The SSCF thermodynamics methods defined here dispatch on
# `ThermodynamicsParameters`, so the parameter set itself serves as the `thermodynamics_backend``

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD

const TD_version = pkgversion(TD)

# --- version adapters ------------------------------------------------------------------------ #

if TD_version ≥ v"0.11.7" # https://github.com/CliMA/Thermodynamics.jl/pull/191/changes
    @inline _virtual_temperature(ps::TD.Parameters.ThermodynamicsParameters, T, q) = TD.virtual_temperature(ps, T, q)
else
    # < 0.11.7 takes an (unused) positional ρ
    @inline _virtual_temperature(ps::TD.Parameters.ThermodynamicsParameters, T, q) = TD.virtual_temperature(ps, T, zero(T), q)
end

# The 0.15.3 restructure (PR #293) renamed `dry_pottemp_given_pressure` and `q_vap_saturation_generic`.
if TD_version < v"0.15.3"
    @inline _dry_pottemp(ps::TD.Parameters.ThermodynamicsParameters, T, p) = TD.dry_pottemp_given_pressure(ps, T, p)
    @inline _q_vap_saturation_ρ(ps::TD.Parameters.ThermodynamicsParameters, T, ρ, phase) = TD.q_vap_saturation_generic(ps, T, ρ, phase)
else
    @inline _dry_pottemp(ps::TD.Parameters.ThermodynamicsParameters, T, p) = TD.potential_temperature_given_pressure(ps, T, p)
    @inline _q_vap_saturation_ρ(ps::TD.Parameters.ThermodynamicsParameters, T, ρ, phase) = TD.q_vap_saturation(ps, T, ρ, phase)
end

@inline _partition(q_tot, q_liq, q_ice) = TD.PhasePartition(q_tot, q_liq, q_ice)
@inline _partition(q_tot) = TD.PhasePartition(q_tot)

# --- physical constants ------------------------------------------------------------------------ #

@inline SSCF.R_d(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params))
@inline SSCF.R_v(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_v(thermo_params))
@inline SSCF.grav(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.grav(thermo_params))
# ε = M_v/M_d ≈ 0.622 = R_d/R_v. Thermodynamics' own `molmass_ratio` is the INVERSE (M_d/M_v ≈ 1.608),
# so compute ε directly from R_d/R_v — matches the backend contract and is version-independent.
@inline SSCF.molmass_ratio(thermo_params::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.R_d(thermo_params) / TD.Parameters.R_v(thermo_params))
@inline SSCF.cp_d(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.cp_d(p))
@inline SSCF.cp_v(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.cp_v(p))
@inline SSCF.cp_l(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.cp_l(p))
@inline SSCF.cp_i(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.cp_i(p))
@inline SSCF.p_ref(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.p_ref_theta(p))
@inline SSCF.T_0(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.T_0(p))
@inline SSCF.T_freeze(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.T_freeze(p))
@inline SSCF.T_icenuc(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.T_icenuc(p))
@inline SSCF.L_v0(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.LH_v0(p))
@inline SSCF.L_s0(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.LH_s0(p))
@inline SSCF.e_ref(p::TD.Parameters.ThermodynamicsParameters{FT}, ::Type{FT2} = FT) where {FT, FT2} = FT2(TD.Parameters.press_triple(p))

@inline SSCF.latent_heat_generic(thermo_params::TD.Parameters.ThermodynamicsParameters, T, LH_0, Δcp) = TD.latent_heat_generic(thermo_params, T, LH_0, Δcp)
@inline SSCF.latent_heat_vapor(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = TD.latent_heat_vapor(thermo_params, T)
@inline SSCF.latent_heat_sublim(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = TD.latent_heat_sublim(thermo_params, T)

# --- saturation vapor pressure ------------------------------------------------------------------ #

@inline SSCF.saturation_vapor_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ::SSCF.Liquid) = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
@inline SSCF.saturation_vapor_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ::SSCF.Ice) = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
@inline SSCF.saturation_vapor_pressure_liq(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = SSCF.saturation_vapor_pressure(thermo_params, T, SSCF.Liquid())
@inline SSCF.saturation_vapor_pressure_ice(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = SSCF.saturation_vapor_pressure(thermo_params, T, SSCF.Ice())

"""
    _p_v_sat_mixture(thermo_params, T, λ)

Saturation vapor pressure over a liquid/ice mixture at liquid fraction `λ`: the Rankine–Kirchhoff
expression with the reference latent heat and heat-capacity difference weighted by `λ`.
"""
@inline function _p_v_sat_mixture(thermo_params::TD.Parameters.ThermodynamicsParameters, T, λ)
    LH_v0 = TD.Parameters.LH_v0(thermo_params)
    LH_s0 = TD.Parameters.LH_s0(thermo_params)
    cp_v = TD.Parameters.cp_v(thermo_params)
    cp_l = TD.Parameters.cp_l(thermo_params)
    cp_i = TD.Parameters.cp_i(thermo_params)
    LH_0 = λ * LH_v0 + (one(λ) - λ) * LH_s0
    Δcp = λ * (cp_v - cp_l) + (one(λ) - λ) * (cp_v - cp_i)
    return TD.saturation_vapor_pressure(thermo_params, T, LH_0, Δcp)
end

# The temperature ramp between `T_icenuc` and `T_freeze`, which is the branch `PhaseEquil` selects.
@inline SSCF.liquid_fraction(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = TD.liquid_fraction(thermo_params, T, TD.PhaseEquil)

@inline SSCF.q_vap_saturation_from_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, q_tot, p, T) =
    TD.q_vap_saturation_from_pressure(thermo_params, q_tot, p, T, TD.PhaseEquil)

@inline function _q_sat_from_e(thermo_params::TD.Parameters.ThermodynamicsParameters, e_sat, p)
    ε = SSCF.molmass_ratio(thermo_params)
    denom = p - (one(ε) - ε) * e_sat
    return denom > zero(denom) ? ε * e_sat / denom : one(denom)
end

@inline SSCF.q_vap_saturation(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, phase::SSCF.AbstractPhase) =
    _q_sat_from_e(thermo_params, SSCF.saturation_vapor_pressure(thermo_params, T, phase), p)
@inline SSCF.q_vap_saturation(
    thermo_params::TD.Parameters.ThermodynamicsParameters, T, p;
    λ = SSCF.liquid_fraction(thermo_params, T),
) = _q_sat_from_e(thermo_params, _p_v_sat_mixture(thermo_params, T, λ), p)

SSCF.q_vap_saturation_liq(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    _q_vap_saturation_ρ(thermo_params, T, TD.air_density(thermo_params, T, p), TD.Liquid())

SSCF.q_vap_saturation_ice(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    _q_vap_saturation_ρ(thermo_params, T, TD.air_density(thermo_params, T, p), TD.Ice())

# --- equilibrium condensate partition from (T, p, q_tot) ---------------------------------------- #

function SSCF.equilibrium_condensate(
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    T,
    p,
    q_tot;
    λ = SSCF.liquid_fraction(thermo_params, T),
)
    R_v = TD.Parameters.R_v(thermo_params)
    ρ = TD.air_density(thermo_params, T, p, _partition(q_tot))
    q_vap_sat = _p_v_sat_mixture(thermo_params, T, λ) / (ρ * R_v * T)
    q_c = max(zero(q_tot), q_tot - q_vap_sat)
    return (; q_liq = λ * q_c, q_ice = (one(λ) - λ) * q_c)
end

# ------------------------------------------------------------ #

SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice) =
    TD.air_density(thermo_params, T, p, _partition(q_tot, q_liq, q_ice))
# (backend, T, p, q_tot) -> ρ — moist density from total water (the form the pipeline calls).
SSCF.air_density(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot) =
    TD.air_density(thermo_params, T, p, _partition(q_tot))

function SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.virtual_temperature(thermo_params, T, q_tot, q_liq, q_ice)
end

# ------------------------------------------------------------ #

# (backend, T, q_tot, q_liq, q_ice) -> T_v — matches the default backend's signature (no `p`; the
# partition suffices). `lev_to_z` calls this 4-argument (after-backend) form.
SSCF.virtual_temperature(thermo_params::TD.Parameters.ThermodynamicsParameters, T, q_tot, q_liq, q_ice) =
    _virtual_temperature(thermo_params, T, _partition(q_tot, q_liq, q_ice))

# ------------------------------------------------------------ #

# liquid-ice potential temperature, computed from pressure with the resolved condensate partition
# (`liquid_ice_pottemp_given_pressure` takes `p`, not density).
SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, q_liq, q_ice) =
    TD.liquid_ice_pottemp_given_pressure(thermo_params, T, p, _partition(q_tot, q_liq, q_ice))
function SSCF.liquid_ice_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p, q_tot)
    (; q_liq, q_ice) = SSCF.equilibrium_condensate(thermo_params, T, p, q_tot)
    return SSCF.liquid_ice_pottemp(thermo_params, T, p, q_tot, q_liq, q_ice)
end

# ------------------------------------------------------------ #

SSCF.dry_pottemp(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) = _dry_pottemp(thermo_params, T, p)

# ------------------------------------------------------------ #

@inline function q_vap_saturation_fixed_liquid_fraction(param_set::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, λ)
    R_d = TD.Parameters.R_d(param_set)
    R_v = TD.Parameters.R_v(param_set)
    FT = eltype(param_set)
    p_v_sat = _p_v_sat_mixture(param_set, T, λ)
    Δp = p - p_v_sat
    return Δp ≥ eps(FT) ? R_d / R_v * (1 - q_tot) * p_v_sat / Δp : one(Δp) / eps(FT)
end

"""
    _θ_li_derivative_state(param_set, T, p, q_tot, vars)

Quantities needed to differentiate `θ_li` with respect to `T` at fixed pressure, given the
condensate humidities and their temperature derivatives in `vars`.
"""
@inline function _θ_li_derivative_state(param_set::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, vars)
    R_v = TD.Parameters.R_v(param_set)
    cp_v = TD.Parameters.cp_v(param_set)
    cp_l = TD.Parameters.cp_l(param_set)
    cp_i = TD.Parameters.cp_i(param_set)
    LH_v0 = TD.Parameters.LH_v0(param_set)
    LH_s0 = TD.Parameters.LH_s0(param_set)
    p0 = TD.Parameters.p_ref_theta(param_set)

    q = _partition(q_tot, vars.q_liq, vars.q_ice)
    R_m = TD.gas_constant_air(param_set, q)
    cp_m = TD.cp_m(param_set, q)
    α = R_m / cp_m

    ln_p_over_p0 = log(p / p0)
    θ = T / exp(α * ln_p_over_p0)  # Π = (p/p₀)^α

    L_c = LH_v0 * vars.q_liq + LH_s0 * vars.q_ice
    F = 1 - L_c / (cp_m * T)

    ∂R_m_∂T = R_v * vars.∂qvs_∂T
    ∂cp_m_∂T = (cp_l - cp_v) * vars.∂q_liq_∂T + (cp_i - cp_v) * vars.∂q_ice_∂T
    ∂α_∂T = (∂R_m_∂T * cp_m - R_m * ∂cp_m_∂T) / cp_m^2

    ∂L_c_∂T = LH_v0 * vars.∂q_liq_∂T + LH_s0 * vars.∂q_ice_∂T
    ∂F_∂T = -1 / (cp_m * T) * (∂L_c_∂T - L_c * (1 / T + ∂cp_m_∂T / cp_m))

    return (; θ, F, ln_p_over_p0, ∂α_∂T, ∂F_∂T)
end

@inline function _fixed_λ_θ_li_and_derivative(param_set::TD.Parameters.ThermodynamicsParameters, T, p, q_tot, λ)
    q_vap_sat = q_vap_saturation_fixed_liquid_fraction(param_set, T, p, q_tot, λ)

    q_cond = max(q_tot - q_vap_sat, zero(q_tot))
    saturated = q_cond > zero(q_tot)

    q_liq = λ * q_cond
    q_ice = (1 - λ) * q_cond

    if saturated
        FT = eltype(param_set)
        R_v = TD.Parameters.R_v(param_set)
        T_0 = TD.Parameters.T_0(param_set)
        LH_v0 = TD.Parameters.LH_v0(param_set)
        LH_s0 = TD.Parameters.LH_s0(param_set)
        cp_v = TD.Parameters.cp_v(param_set)
        cp_l = TD.Parameters.cp_l(param_set)
        cp_i = TD.Parameters.cp_i(param_set)

        p_v_sat = _p_v_sat_mixture(param_set, T, λ)
        Δp = p - p_v_sat

        # Latent heat of the mixture at `T`, the temperature derivative of the Rankine–Kirchhoff
        # exponent: d(ln p*)/dT = L(T) / (R_v T²).
        LH_0 = λ * LH_v0 + (1 - λ) * LH_s0
        Δcp = λ * (cp_v - cp_l) + (1 - λ) * (cp_v - cp_i)
        ∂lnp_v_sat_∂T = (LH_0 + Δcp * (T - T_0)) / (R_v * T^2)

        amplification = ifelse(Δp ≥ eps(FT), p / Δp, zero(Δp))
        ∂qvs_∂T = q_vap_sat * ∂lnp_v_sat_∂T * amplification

        # λ is prescribed, so ∂λ/∂T = 0.
        ∂q_liq_∂T = -λ * ∂qvs_∂T
        ∂q_ice_∂T = -(1 - λ) * ∂qvs_∂T
    else
        # q_liq = q_ice = 0 because the max() has clamped the condensate to zero, so the state has
        # no dependence on q_vap_sat.
        ∂qvs_∂T = zero(q_tot)
        ∂q_liq_∂T = zero(q_tot)
        ∂q_ice_∂T = zero(q_tot)
    end

    vars = (; q_liq, q_ice, ∂qvs_∂T, ∂q_liq_∂T, ∂q_ice_∂T)

    st = _θ_li_derivative_state(param_set, T, p, q_tot, vars)

    θ_li_val = TD.liquid_ice_pottemp_given_pressure(param_set, T, p, _partition(q_tot, q_liq, q_ice))

    # Fixed-pressure product rule for θ_li = θ·F, with λ held fixed.
    ∂θ_∂T = st.θ * (1 / T - st.ln_p_over_p0 * st.∂α_∂T)
    ∂θ_li_∂T = ∂θ_∂T * st.F + st.θ * st.∂F_∂T

    return θ_li_val, ∂θ_li_∂T, q_liq, q_ice
end

# Temperature of unsaturated air at pressure `p` with liquid-ice potential temperature `θ_li`:
# with no condensate θ_li reduces to the potential temperature of the moist mixture.
@inline function _unsaturated_temperature(param_set::TD.Parameters.ThermodynamicsParameters, p, θ_li, q_tot)
    q = _partition(q_tot)
    α = TD.gas_constant_air(param_set, q) / TD.cp_m(param_set, q)
    return θ_li * (p / TD.Parameters.p_ref_theta(param_set))^α
end

function saturation_adjustment_given_liquid_fraction(::Type{TD.RS.NewtonsMethod}, param_set::TD.Parameters.ThermodynamicsParameters, p, θ_li, q_tot, λ, maxiter, tol)
    FT = eltype(param_set)

    T_unsat = _unsaturated_temperature(param_set, p, θ_li, q_tot)

    # Check saturation using the prescribed λ.
    q_vap_unsat = q_vap_saturation_fixed_liquid_fraction(param_set, T_unsat, p, q_tot, λ)

    if q_tot <= q_vap_unsat
        return (; T = T_unsat, q_liq = zero(q_tot), q_ice = zero(q_tot), converged = true)
    end

    T_floor = sqrt(eps(FT))
    T_guess = max(T_unsat, T_floor)

    roots_function = T -> begin
        T_val = max(T, T_floor)
        θ_li_val, ∂θ_li_∂T, _, _ = _fixed_λ_θ_li_and_derivative(param_set, T_val, p, q_tot, λ)
        (θ_li_val - θ_li, ∂θ_li_∂T)
    end

    sol = TD.RS.find_zero(
        roots_function,
        TD.RS.NewtonsMethod(T_guess),
        TD.solution_type(),
        tol isa TD.RS.AbstractTolerance ? tol : TD.RS.RelativeSolutionTolerance(tol),
        maxiter,
    )
    T, converged = sol.root, sol.converged

    _, _, q_liq, q_ice = _fixed_λ_θ_li_and_derivative(param_set, T, p, q_tot, λ)

    return (; T, q_liq, q_ice, converged)
end

# `λ = NaN` (the default backend's sentinel) means "let Thermodynamics partition by temperature";
# any other value holds the liquid fraction fixed.
function SSCF.saturation_adjust_pθq(
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    p,
    θ_liq_ice,
    q_tot;
    maxiter::Int = 50,
    tol = nothing,
    λ = NaN,
)
    FT = eltype(thermo_params)
    if isnan(λ)
        # `PhaseEquil_pθq` runs the θ_liq_ice saturation adjustment internally. `tol === nothing`
        # lets Thermodynamics use its own default relative tolerance.
        _tol = tol === nothing ? nothing : FT(tol)
        ts = TD.PhaseEquil_pθq(thermo_params, FT(p), FT(θ_liq_ice), FT(q_tot), maxiter, _tol)
        q = TD.PhasePartition(thermo_params, ts)
        return (; T = TD.air_temperature(thermo_params, ts), q_liq = q.liq, q_ice = q.ice)
    end
    _tol = tol === nothing ? FT(1.0e-6) : FT(tol)
    sat = saturation_adjustment_given_liquid_fraction(
        TD.RS.NewtonsMethod, thermo_params, FT(p), FT(θ_liq_ice), FT(q_tot), FT(λ), maxiter, _tol,
    )
    return (; T = sat.T, q_liq = sat.q_liq, q_ice = sat.q_ice)
end

# --- surface quantities -------------------------------------------------------------------------- #

function SSCF.saturation_specific_humidity_from_pT(thermo_params::TD.Parameters.ThermodynamicsParameters, p, T, phase::TD.Phase = TD.Liquid())
    ρg = TD.air_density(thermo_params, T, p)
    return _q_vap_saturation_ρ(thermo_params, T, ρg, phase)  # surface specific humidity over liquid
end

function SSCF.saturation_mixing_ratio_from_pT(thermo_params::TD.Parameters.ThermodynamicsParameters, p, T, phase::TD.Phase = TD.Liquid())
    pv = TD.saturation_vapor_pressure(thermo_params, T, phase)
    ε = SSCF.molmass_ratio(thermo_params)  # M_v/M_d ≈ 0.622 (see the accessor above; NOT TD's inverse)
    return ε * pv / (p - pv)  # saturation total-water mixing ratio at the surface: w_s = ε·e_s/(p−e_s)
end

end # module
