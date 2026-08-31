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
@inline SSCF.latent_heat_vapor(thermo_params::TD.Parameters.ThermodynamicsParameters, T) =  TD.latent_heat_vapor(thermo_params, T)
@inline SSCF.latent_heat_sublim(thermo_params::TD.Parameters.ThermodynamicsParameters, T) =  TD.latent_heat_sublim(thermo_params, T)


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

@inline SSCF.saturation_vapor_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ::SSCF.Liquid) = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
@inline SSCF.saturation_vapor_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, T, ::SSCF.Ice) = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
@inline SSCF.saturation_vapor_pressure_liq(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = SSCF.saturation_vapor_pressure(thermo_params, T, SSCF.Liquid())
@inline SSCF.saturation_vapor_pressure_ice(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = SSCF.saturation_vapor_pressure(thermo_params, T, SSCF.Ice())

@inline SSCF.liquid_fraction(thermo_params::TD.Parameters.ThermodynamicsParameters, T) = TD.liquid_fraction(thermo_params, T, zero(T), zero(T))

@inline SSCF.q_vap_saturation_from_pressure(thermo_params::TD.Parameters.ThermodynamicsParameters, q_tot, p, T) = TD.q_vap_saturation_from_pressure(thermo_params, q_tot, p, T)

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
) = _q_sat_from_e(thermo_params, TD.saturation_vapor_pressure_mixture(thermo_params, T, λ), p)

# ------------------------------------------------------------ #


@inline function q_vap_saturation_fixed_liquid_fraction(
    param_set::TD.Parameters.ThermodynamicsParameters,
    T,
    p,
    q_tot,
    λ,
)
    p_v_sat = TD.saturation_vapor_pressure_mixture(param_set, T, λ)
    return TD.q_vap_saturation_from_pressure_calc(
        param_set,
        q_tot,
        p,
        p_v_sat,
    )
end

@inline function _fixed_λ_θ_li_and_derivative(
    param_set::TD.Parameters.ThermodynamicsParameters,
    T,
    p,
    q_tot,
    λ,
)
    q_vap_sat = q_vap_saturation_fixed_liquid_fraction(
        param_set,
        T,
        p,
        q_tot,
        λ,
    )

    q_cond = max(q_tot - q_vap_sat, zero(q_tot))
    saturated = q_cond > zero(q_tot)

    q_liq = λ * q_cond
    q_ice = (1 - λ) * q_cond

    if saturated
        FT = eltype(param_set)
        R_v = TD.TP.R_v(param_set)

        p_v_sat = TD.saturation_vapor_pressure_mixture(param_set, T, λ)

        Δp = p - p_v_sat
        ∂lnp_v_sat_∂T =
            TD.latent_heat_mixed(param_set, T, λ) / (R_v * T^2)

        amplification = ifelse(
            Δp ≥ TD.ϵ_numerics(FT),
            p / Δp,
            zero(Δp),
        )

        ∂qvs_∂T =
            q_vap_sat * ∂lnp_v_sat_∂T * amplification

        # λ is prescribed, so ∂λ/∂T = 0.
        ∂q_liq_∂T = -λ * ∂qvs_∂T
        ∂q_ice_∂T = -(1 - λ) * ∂qvs_∂T
    else
        # q_liq = q_ice = 0 because the max() has clamped the
        # condensate to zero. Therefore the state itself has no
        # dependence on q_vap_sat.
        ∂qvs_∂T = zero(q_tot)
        ∂q_liq_∂T = zero(q_tot)
        ∂q_ice_∂T = zero(q_tot)
    end

    vars = (;
        q_liq,
        q_ice,
        ∂qvs_∂T,
        ∂q_liq_∂T,
        ∂q_ice_∂T,
    )

    st = TD._θ_li_derivative_state(
        param_set,
        T,
        p,
        q_tot,
        vars,
    )

    θ_li_val = TD.liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        p,
        q_tot,
        q_liq,
        q_ice,
    )

    # This is the same fixed-pressure product rule used by
    # Thermodynamics.jl's ∂θ_li_∂T_sat_p, except λ is fixed here.
    ∂θ_∂T =
        st.θ * (
            1 / T -
            st.ln_p_over_p0 * st.∂α_∂T
        )

    ∂θ_li_∂T =
        ∂θ_∂T * st.F +
        st.θ * st.∂F_∂T

    return θ_li_val, ∂θ_li_∂T, q_liq, q_ice
end



function saturation_adjustment_given_liquid_fraction(::Type{TD.RS.NewtonsMethod}, param_set::TD.Parameters.ThermodynamicsParameters, ::TD.pθ_li, p, θ_li, q_tot, λ, maxiter, tol)
    FT = eltype(param_set)

    # Temperature with no condensate.
    T_unsat = TD.air_temperature(param_set, TD.pθ_li(), p, θ_li, q_tot)

    # Check saturation using the prescribed λ.
    q_vap_unsat = q_vap_saturation_fixed_liquid_fraction(param_set, T_unsat, p, q_tot, λ)

    if q_tot <= q_vap_unsat
        return (; T = T_unsat, q_liq = zero(q_tot), q_ice = zero(q_tot), converged = true)
    end

    T_guess = max(T_unsat, TD.T_positive_floor(FT))

    roots_function = T -> begin
        T_val = max(T, TD.T_positive_floor(FT))
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
    _tol = tol === nothing ? FT(1.0e-6) : FT(tol)
    sat = if isnan(λ)
        TD.saturation_adjustment(TD.RS.NewtonsMethod, thermo_params, TD.pθ_li(), FT(p), FT(θ_liq_ice), FT(q_tot), maxiter, _tol)
    else
        saturation_adjustment_given_liquid_fraction(
            TD.RS.NewtonsMethod, thermo_params, TD.pθ_li(), FT(p), FT(θ_liq_ice), FT(q_tot), FT(λ), maxiter, _tol,
        )
    end
    return (; T = sat.T, q_liq = sat.q_liq, q_ice = sat.q_ice)
end

SSCF.q_vap_saturation_liq(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    TD.q_vap_saturation(thermo_params, T, TD.air_density(thermo_params, T, p), TD.Liquid())

SSCF.q_vap_saturation_ice(thermo_params::TD.Parameters.ThermodynamicsParameters, T, p) =
    TD.q_vap_saturation(thermo_params, T, TD.air_density(thermo_params, T, p), TD.Ice())

function SSCF.saturation_specific_humidity_from_pT(thermo_params::TD.Parameters.ThermodynamicsParameters, p, T, phase::TD.Phase = TD.Liquid())
    ρg = TD.air_density(thermo_params, T, p)
    return TD.q_vap_saturation(thermo_params, T, ρg, phase)  # surface specific humidity over liquid
end


function SSCF.saturation_mixing_ratio_from_pT(thermo_params::TD.Parameters.ThermodynamicsParameters, p, T, phase::TD.Phase = TD.Liquid())
    pv = TD.saturation_vapor_pressure(thermo_params, T, phase)
    ε = SSCF.molmass_ratio(thermo_params)  # M_v/M_d ≈ 0.622 (see the accessor above; NOT TD's inverse)
    return ε * pv / (p - pv)  # saturation total-water mixing ratio at the surface: w_s = ε·e_s/(p−e_s)
end



""" Get Thermodynamics.jl parameter set from a NamedTuple of parameters (default value if not passed). """
function get_thermo_params(params::NamedTuple)
    return TD.Parameters.ThermodynamicsParameters(
        error("not impleemnted yet")
    )
end


end # module
