#=
Physical diagnostics computed from the Atlas archive rather than read from it. The registry in
`units.jl` calls these from its derived specs; they take SI arguments and return SI results.
=#

"""Atlas's `QSMALL`: the specific humidity below which a species counts as absent."""
const ATLAS_Q_SMALL = 1.0e-14

"""
Latent heat of vaporization [J/kg]
"""
const ATLAS_LH_V0 = 2_510_400.0

"""
Latent heat of sublimation [J/kg]
"""
const ATLAS_LH_S0 = 2_844_000.0

"""
Dry-air heat capacity [J/kg/K]
"""
const ATLAS_CP_D = 1004.0

"""Terminal-velocity caps [m/s] for velocities diagnosed from a sedimentation flux."""
const ATLAS_SED_W_CAP = (; cloud = -0.2, ice = -1.0, snow = -2.0, rain = -10.0)

"""Threshold for cloud-ice autoconversion as a radius [m]; Morrison splits at diameter `DCS = 125 µm`."""
const ATLAS_ICE_THRESHOLD_RADIUS = 62.5e-6

"""Smallest particle radius resolved [m]; regularizes the size distribution slope as `q` goes to zero."""
const ATLAS_PARTICLE_MIN_RADIUS = 0.2e-6

"""
Apparent density of cloud ice [kg/m3] in the M2005 microphysics the archive was written with.

This is the apparent density of the ice particles the size distribution describes, not the bulk
density of solid ice (917 kg/m3).
"""
const ATLAS_ICE_APPARENT_DENSITY = 500.0

"""
    sedimentation_velocity(F, q, ρ, cap)

Mass-weighted fall speed [m/s] implied by a sedimentation flux `F` [kg/m2/s] of a species with
specific content `q` and air density `ρ`, positive downward.

Where the species is absent `F/(qρ)` is `0/0` and is reported as zero fall speed; where `q` is tiny
but non-zero the quotient is unbounded, so it is clamped into `[cap, 0]`.
"""
function sedimentation_velocity(F::FT, q::FT, ρ::FT, cap::FT) where {FT}
    w = F / (q * ρ)
    isfinite(w) || return zero(FT)
    return -clamp(w, cap, zero(FT))
end

"""
    ice_distribution_slope(q, N, ρ, ρ_ice, μ = 0; r_min = ATLAS_PARTICLE_MIN_RADIUS)

Slope `λ` [1/m] of the gamma ice size distribution implied by specific content `q`, number per
volume `N` and air density `ρ`, with shape parameter `μ`.

For `n(r) = n₀ r^μ e^{-λr}` and `m(r) = (4/3)π ρ_ice r³`, `N = n₀ Γ(μ+1)/λ^{μ+1}` and
`ρq = (4/3)π ρ_ice n₀ Γ(μ+4)/λ^{μ+4}`, so `λ³ = (4/3)π ρ_ice (μ+1)(μ+2)(μ+3) N/(ρq)` using
`Γ(μ+4)/Γ(μ+1) = (μ+1)(μ+2)(μ+3)`, exact for any real `μ`. Returns `NaN` where there is no ice.
"""
function ice_distribution_slope(
    q::FT,
    N::FT,
    ρ::FT,
    ρ_ice::FT,
    μ::FT = zero(FT);
    r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS),
) where {FT}
    (isfinite(q) && isfinite(N) && isfinite(ρ) && q > 0 && N > 0) || return FT(NaN)
    m_min = FT(4 // 3) * FT(π) * ρ_ice * r_min^3
    q_eff = q + N * m_min / ρ
    shape = (μ + one(FT)) * (μ + FT(2)) * (μ + FT(3))
    return cbrt(FT(4 // 3) * FT(π) * ρ_ice * shape * N / (ρ * q_eff))
end

"""
    ice_mean_radius(q, N, ρ, ρ_ice, μ = 0; r_min = ATLAS_PARTICLE_MIN_RADIUS)

Number-weighted mean ice radius [m], `⟨r⟩ = ∫r n dr / ∫n dr = (μ+1)/λ`.
"""
ice_mean_radius(
    q::FT,
    N::FT,
    ρ::FT,
    ρ_ice::FT,
    μ::FT = zero(FT);
    r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS),
) where {FT} = (μ + one(FT)) / ice_distribution_slope(q, N, ρ, ρ_ice, μ; r_min)

"""
    ice_process_threshold_weight(q, N, ρ, ρ_ice; r_threshold = ATLAS_ICE_THRESHOLD_RADIUS)

Deposition-weighted fraction of the ice distribution below `r_threshold`, `1 - e^{-λr}(1 + λr)`.

The radius weighting is the one deposition takes, since `dm/dt = S/τ` with `τ ∝ r`. This closed form
holds for the exponential distribution, which is the Atlas LES ice distribution (`N0I = NI·λᵢ`), so
`μ` is fixed at zero rather than exposed.
"""
function ice_process_threshold_weight(
    q::FT,
    N::FT,
    ρ::FT,
    ρ_ice::FT;
    r_threshold::FT = FT(ATLAS_ICE_THRESHOLD_RADIUS),
    kwargs...,
) where {FT}
    λ = ice_distribution_slope(q, N, ρ, ρ_ice; kwargs...)
    x = λ * r_threshold
    return one(FT) - exp(-x) * (one(FT) + x)
end

"""
    water_vapor_diffusivity_in_air(T, p)

Water-vapour diffusivity in air [m2/s] at temperature `T` [K] and pressure `p` [Pa]
(Pruppacher and Klett, 1997).
"""
water_vapor_diffusivity_in_air(T::FT, p::FT) where {FT} =
    FT(2.11e-5) * (T / FT(273.15))^FT(1.94) * (FT(101325) / p)

"""
    ice_theory_timescale(q, N, ρ, T, p, ρ_ice, μ = 0; r_min = ATLAS_PARTICLE_MIN_RADIUS)

Supersaturation relaxation timescale for ice [s], `1/(4π D N r ρ)`, with the mean radius implied by
`(q, N, ρ)` and the vapour diffusivity at `(T, p)`.
"""
function ice_theory_timescale(
    q::FT,
    N::FT,
    ρ::FT,
    T::FT,
    p::FT,
    ρ_ice::FT,
    μ::FT = zero(FT);
    r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS),
) where {FT}
    r = ice_mean_radius(q, N, ρ, ρ_ice, μ; r_min)
    D = water_vapor_diffusivity_in_air(T, p)
    return inv(4 * FT(π) * D * N * r * ρ)
end

"""
    cooper_ice_nucleus_concentration(T)

Cooper (1986) ice-nucleus concentration [1/m3] at temperature `T` [K].
"""
cooper_ice_nucleus_concentration(T::FT) where {FT} =
    FT(0.005) * exp(FT(0.304) * (FT(273.15) - T)) * FT(1000)
