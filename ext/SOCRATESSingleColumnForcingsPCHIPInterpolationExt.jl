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

function SSCF.Interpolation.conservative_mass_matrix(::PCHIPInterpolationMethod)
    error("Not implemented")
end

function SSCF.Interpolation.conservative_mass_matrix(::PCHIPSmoothDerivativeInterpolationMethod)
    error("Not implemented")
end

function SSCF.Interpolation.build_spline(
    ::PCHIPInterpolationMethod,
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
)
    spl = PCHIPInterpolation.Interpolator(xp, fp)
    if bc isa SSCF.Interpolation.ErrorBoundaryCondition
        error("Not implemented: pchip with bc = error")
    elseif bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
        return pchip_extrapolate(spl)
    else
        error("Unsupported boundary condition for PCHIP: $bc")
    end
end

function SSCF.Interpolation.build_spline(
    method::PCHIPSmoothDerivativeInterpolationMethod,
    xp,
    fp;
    bc::SSCF.Interpolation.ValidBoundaryConditions = SSCF.Interpolation.ErrorBoundaryCondition(),
)
    return pchip_smooth_derivative(
        xp,
        fp;
        bc = bc,
        f_enhancement_factor = method.f_enhancement_factor,
        f_p_enhancement_factor = method.f_p_enhancement_factor,
    )
end

# Continue the edge slope beyond the data for a continuous (not necessarily smooth) derivative.
function pchip_extrapolate(spl::PCHIPInterpolation.Interpolator)
    xmin = spl.xs[1]
    xmax = spl.xs[end]
    ymin = spl.ys[1]
    ymax = spl.ys[end]
    return x -> begin
        if x < xmin
            dydx_xmin = PCHIPInterpolation._derivative(spl, Val(:begin), 1)
            return ymin - (xmin - x) * dydx_xmin
        elseif x > xmax
            dydx_xmax = PCHIPInterpolation._derivative(spl, Val(:end), length(spl.xs) - 1)
            return ymax + (x - xmax) * dydx_xmax
        else
            return spl(x)
        end
    end
end

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

end # module
