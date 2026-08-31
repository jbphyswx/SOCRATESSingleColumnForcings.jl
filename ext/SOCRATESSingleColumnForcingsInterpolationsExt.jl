module SOCRATESSingleColumnForcingsInterpolationsExt

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Interpolations: Interpolations

# parametric over the Interpolations.jl spec itself, e.g.
# `InterpolationsInterpolationMethod(Interpolations.Gridded(Interpolations.Linear()))`
struct InterpolationsInterpolationMethod{D} <: SSCF.Interpolation.AbstractInterpolationMethod
    degree::D
end


function SSCF.Interpolation.build_spline(
    method::InterpolationsInterpolationMethod,
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
    drop_collinear::Val = Val(false),
    collinear_tol = promote_type(eltype(xp), eltype(fp))(NaN),
)
    xp, fp = SSCF.Interpolation._maybe_prune(drop_collinear, xp, fp, collinear_tol)
    itp = Interpolations.interpolate((xp,), fp, method.degree)
    extrap =
        bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition ? Interpolations.Line() :
        bc isa SSCF.Interpolation.NearestBoundaryCondition ? Interpolations.Flat() : Interpolations.Throw()
    return Interpolations.extrapolate(itp, extrap)
end

SSCF.Interpolation.interpolant_nodes(etp::Interpolations.AbstractExtrapolation) =
    collect(Interpolations.knots(etp))

SSCF.Interpolation.rebuild_interpolant(etp::Interpolations.AbstractExtrapolation, xs, ys) =
    Interpolations.extrapolate(Interpolations.interpolate((xs,), ys, Interpolations.itpflag(etp)), etp.et)

# Polynomial degree of an Interpolations spec, read off the `Degree{N}` type parameter. A spec that
# is not a `Degree` — a Lanczos kernel, a monotonic interpolator — has no polynomial degree, and no
# fixed-order quadrature is exact for it.
_spline_degree(::Interpolations.Degree{N}) where {N} = N
_spline_degree(spec::Interpolations.BSpline) = _spline_degree(spec.degree)
_spline_degree(spec::Interpolations.Gridded) = _spline_degree(spec.degree)
_spline_degree(spec) = error(
    "safe_integrate: the Interpolations spec $(spec) carries no polynomial degree, so no " *
    "fixed-order quadrature integrates it exactly. Use a BSpline or Gridded spec built from " *
    "Constant, Linear, Quadratic or Cubic.",
)

# Gauss-Legendre nodes and weights on [-1, 1]; an n-point rule is exact to degree 2n - 1.
_gauss_rule(::Val{1}) = ((0.0,), (2.0,))
_gauss_rule(::Val{2}) = ((-0.5773502691896257, 0.5773502691896257), (1.0, 1.0))
_gauss_rule(::Val{3}) = ((-0.7745966692414834, 0.0, 0.7745966692414834), (5 / 9, 8 / 9, 5 / 9))
_gauss_rule(::Val{n}) where {n} = error(
    "safe_integrate: no tabulated Gauss-Legendre rule with $n points (spec degree too high)",
)

"""
    safe_integrate(spl::Interpolations.AbstractExtrapolation, x1, x2; bc)

Integrate an Interpolations-backed spline from `x1` to `x2`.

The spline is a polynomial only BETWEEN knots, so the interval is split at every knot it contains and
a Gauss-Legendre rule of just enough order for the spec's degree is applied to each piece, making the
result exact. The out-of-range policy is baked into the object by `build_spline`, so `bc` is accepted
for interface compatibility and not re-applied.
"""
function SSCF.Interpolation.safe_integrate(
    spl::Interpolations.AbstractExtrapolation,
    x1::FT,
    x2::FT;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
) where {FT}
    x1 == x2 && return zero(FT)
    lo, hi = minmax(x1, x2)
    degree = _spline_degree(Interpolations.itpflag(spl))
    nodes, weights = _gauss_rule(Val(cld(degree + 1, 2)))

    breaks = FT[lo; collect(FT, Interpolations.knotsbetween(spl; start = lo, stop = hi)); hi]
    total = zero(FT)
    @inbounds for i in 1:(length(breaks) - 1)
        a, b = breaks[i], breaks[i + 1]
        b > a || continue
        half, mid = (b - a) / 2, (a + b) / 2
        piece = zero(FT)
        for (ξ, w) in zip(nodes, weights)
            piece += FT(w) * spl(mid + half * FT(ξ))
        end
        total += half * piece
    end
    return x2 < x1 ? -total : total
end

end # module
