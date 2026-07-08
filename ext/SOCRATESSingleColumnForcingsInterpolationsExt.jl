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
)
    itp = Interpolations.interpolate((xp,), fp, method.degree)
    extrap =
        bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition ? Interpolations.Line() :
        bc isa SSCF.Interpolation.NearestBoundaryCondition ? Interpolations.Flat() : Interpolations.Throw()
    return Interpolations.extrapolate(itp, extrap)
end

end # module
