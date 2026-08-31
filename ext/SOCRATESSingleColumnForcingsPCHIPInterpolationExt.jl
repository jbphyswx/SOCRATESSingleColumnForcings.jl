module SOCRATESSingleColumnForcingsPCHIPInterpolationExt

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using PCHIPInterpolation: PCHIPInterpolation
using ForwardDiff: ForwardDiff

abstract type AbstractPCHIPInterpolationMethod <: SSCF.Interpolation.AbstractInterpolationMethod end
struct PCHIPInterpolationMethod <: AbstractPCHIPInterpolationMethod end
Base.@kwdef struct PCHIPSmoothDerivativeInterpolationMethod <: AbstractPCHIPInterpolationMethod
    f_enhancement_factor::Int = 1
    f_p_enhancement_factor::Int = 1
end

function SSCF.Interpolation.build_spline(
    ::PCHIPInterpolationMethod,
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
    drop_collinear::Val = Val(false),
    collinear_tol = promote_type(eltype(xp), eltype(fp))(NaN),
)
    xp, fp = SSCF.Interpolation._maybe_prune(drop_collinear, xp, fp, collinear_tol)
    return PCHIPInterpolant(PCHIPInterpolation.Interpolator(xp, fp), bc)
end

function SSCF.Interpolation.build_spline(
    method::PCHIPSmoothDerivativeInterpolationMethod,
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
    drop_collinear::Val = Val(false),
    collinear_tol = promote_type(eltype(xp), eltype(fp))(NaN),
)
    xp, fp = SSCF.Interpolation._maybe_prune(drop_collinear, xp, fp, collinear_tol)
    return pchip_smooth_derivative(
        xp,
        fp;
        bc = bc,
        f_enhancement_factor = method.f_enhancement_factor,
        f_p_enhancement_factor = method.f_p_enhancement_factor,
    )
end

"""
    PCHIPInterpolant(spl, bc)

PCHIP interpolant carrying its own out-of-range policy: `Extrapolate` continues the edge slope
linearly (a continuous, not necessarily smooth, derivative), `Nearest` holds the edge value, `Error`
throws.
"""
struct PCHIPInterpolant{ITP, T, BCT <: SSCF.Interpolation.ValidBoundaryConditions} <:
       SSCF.Interpolation.AbstractInterpolant
    spl::ITP
    bc::BCT
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    dydx_min::T
    dydx_max::T
end

function PCHIPInterpolant(
    spl::PCHIPInterpolation.Interpolator,
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
)
    dydx_min = PCHIPInterpolation._derivative(spl, Val(:begin), 1)
    dydx_max = PCHIPInterpolation._derivative(spl, Val(:end), length(spl.xs) - 1)
    vals = promote(spl.xs[1], spl.xs[end], spl.ys[1], spl.ys[end], dydx_min, dydx_max)
    return PCHIPInterpolant{typeof(spl), eltype(vals), typeof(bc)}(spl, bc, vals...)
end

# Value of the boundary continuation at `x`, anchored at `(x0, y0)` with edge slope `slope` — the
# pointwise counterpart of `_edge_integral` below.
@inline _edge_value(::SSCF.Interpolation.ExtrapolateBoundaryCondition, y0, slope, x0, x) = y0 + (x - x0) * slope
@inline _edge_value(::SSCF.Interpolation.NearestBoundaryCondition, y0, slope, x0, x) = y0
@inline _edge_value(bc::SSCF.Interpolation.ErrorBoundaryCondition, y0, slope, x0, x) =
    error("requested x = $x lies outside the PCHIP knot range and bc = $bc")

@inline function (s::PCHIPInterpolant)(x)
    x < s.xmin && return _edge_value(s.bc, s.ymin, s.dydx_min, s.xmin, x)
    x > s.xmax && return _edge_value(s.bc, s.ymax, s.dydx_max, s.xmax, x)
    return s.spl(x)
end

SSCF.Interpolation.interpolant_nodes(itp::PCHIPInterpolant) = itp.spl.xs
SSCF.Interpolation.rebuild_interpolant(itp::PCHIPInterpolant, xs, ys) =
    PCHIPInterpolant(PCHIPInterpolation.Interpolator(xs, ys), itp.bc)

"""
Build a function with a smooth second derivative by fitting a pchip spline to the data's first
derivative and integrating it back. This rounds corners rather than hitting the data exactly, but
is smooth. `f_enhancement_factor` densifies the data (linear) before differentiating;
`f_p_enhancement_factor` densifies the derivative before its pchip fit — both reduce rounding.
"""
function pchip_smooth_derivative(
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
)
    if !isone(f_enhancement_factor)
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_enhancement_factor - 1))
        for i in 1:(length(xp) - 1)
            xp_new[((i - 1) * f_enhancement_factor + 1):(i * f_enhancement_factor)] .=
                range(xp[i], stop = xp[i + 1], length = f_enhancement_factor + 1)[1:(end - 1)]
        end
        xp_new[end] = xp[end]
        xp, fp = xp_new, SSCF.Interpolation.interpolate_1d(xp_new, xp, fp, SSCF.Interpolation.FastLinear1DInterpolation; bc = bc)
    end

    f_pchip_spl = PCHIPInterpolation.Interpolator(xp, fp)
    dfdx = ForwardDiff.derivative.(Ref(f_pchip_spl), xp)

    if !isone(f_p_enhancement_factor)
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_p_enhancement_factor - 1))
        for i in 1:(length(xp) - 1)
            xp_new[((i - 1) * f_p_enhancement_factor + 1):(i * f_p_enhancement_factor)] .=
                range(xp[i], stop = xp[i + 1], length = f_p_enhancement_factor + 1)[1:(end - 1)]
        end
        xp_new[end] = xp[end]
        xp, dfdx = xp_new, ForwardDiff.derivative.(Ref(f_pchip_spl), xp_new)
    end

    spl_dfdx = PCHIPInterpolation.Interpolator(xp, dfdx)
    dfp_dx_xmin = PCHIPInterpolation._derivative(spl_dfdx, Val(:begin), 1)
    dfp_dx_xmax = PCHIPInterpolation._derivative(spl_dfdx, Val(:end), length(xp) - 1)
    dfdx_min = spl_dfdx.ys[1]
    dfdx_max = spl_dfdx.ys[end]
    ymin = fp[1]
    xmin = xp[1]

    return x -> begin
        if x < xp[1]
            if bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
                x_0 = xp[1]
                Δx = x_0 - x
                return ymin - (dfdx_min * Δx - dfp_dx_xmin * (+x_0 * Δx - (x_0^2 - x^2) / 2))
            else
                error("Requested x below the spline minimum with bc = $bc; use extrapolate")
            end
        elseif x > xp[end]
            if bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
                x_0 = xp[end]
                Δx = x - x_0
                ymax = ymin + PCHIPInterpolation.integrate(spl_dfdx, xp[1], xp[end])
                return ymax + (dfdx_max * Δx + dfp_dx_xmax * ((x^2 - x_0^2) / 2 - x_0 * Δx))
            else
                error("Requested x above the spline maximum with bc = $bc; use extrapolate")
            end
        else
            return ymin + PCHIPInterpolation.integrate(spl_dfdx, xmin, x)
        end
    end
end

# ∫ over [a, b] of the boundary continuation anchored at `(x0, y0)`: the edge value held constant
# under a nearest bc, or the linear continuation `y0 + slope*(x - x0)` under an extrapolate bc.
_edge_integral(::SSCF.Interpolation.NearestBoundaryCondition, y0, slope, x0, a, b) = y0 * (b - a)
_edge_integral(::SSCF.Interpolation.ExtrapolateBoundaryCondition, y0, slope, x0, a, b) =
    y0 * (b - a) + slope / 2 * ((b - x0)^2 - (a - x0)^2)
_edge_integral(bc::SSCF.Interpolation.ErrorBoundaryCondition, y0, slope, x0, a, b) =
    error("integration range [$a, $b] lies outside the knots and bc = $bc")

"""
    safe_integrate(spl::PCHIPInterpolant, x1, x2; bc = spl.bc)

Integrate the PCHIP interpolant over `[x1, x2]`, honoring `bc` outside the knot range: the interior
uses `PCHIPInterpolation.integrate`, and each out-of-range part is integrated in closed form from the
edge value and edge slope. Reversed bounds give the negated integral.
"""
function SSCF.Interpolation.safe_integrate(
    s::PCHIPInterpolant,
    x1::T1,
    x2::T2;
    bc::SSCF.Interpolation.ValidBoundaryConditions = s.bc,
) where {T1 <: Real, T2 <: Real}
    FT = promote_type(T1, T2, typeof(s.xmin), typeof(s.ymin))
    lo, hi = FT(min(x1, x2)), FT(max(x1, x2))
    if (bc isa SSCF.Interpolation.ErrorBoundaryCondition) && ((lo < s.xmin) || (hi > s.xmax))
        error("integration bounds [$lo, $hi] fall outside the PCHIP knot range [$(s.xmin), $(s.xmax)] and bc = $bc")
    end

    y = zero(FT)
    lo < s.xmin && (y += _edge_integral(bc, s.ymin, s.dydx_min, s.xmin, lo, min(hi, s.xmin)))
    a, b = max(lo, s.xmin), min(hi, s.xmax)
    b > a && (y += PCHIPInterpolation.integrate(s.spl, a, b))
    hi > s.xmax && (y += _edge_integral(bc, s.ymax, s.dydx_max, s.xmax, max(lo, s.xmax), hi))

    return (x2 < x1) ? -y : y
end

end # module
