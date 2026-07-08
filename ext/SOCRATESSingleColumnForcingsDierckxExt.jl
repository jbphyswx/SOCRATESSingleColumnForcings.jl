module SOCRATESSingleColumnForcingsDierckxExt

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Dierckx: Dierckx

struct DierckxSpline1DInterpolationMethod <: SSCF.Interpolation.AbstractInterpolationMethod
    k::Int
end

function SSCF.Interpolation.build_spline(
    method::DierckxSpline1DInterpolationMethod,
    xp,
    fp;
    bc::BCT = SSCF.Interpolation.ErrorBoundaryCondition(),
) where {BCT <: SSCF.Interpolation.ValidBoundaryConditions}
    return Dierckx.Spline1D(xp, fp; k = method.k, bc = SSCF.Interpolation.bc_string(bc))
end

function SSCF.Interpolation.interpolate_1d(
    x,
    xp,
    fp,
    method::DierckxSpline1DInterpolationMethod;
    bc::BCT = SSCF.Interpolation.ErrorBoundaryCondition(),
) where {BCT <: SSCF.Interpolation.ValidBoundaryConditions}
    spl = SSCF.Interpolation.build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end

"""
Dierckx does not respect boundary-condition mode (extrapolate/nearest/error) when integrating,
so integration outside the knots is handled here (unlike Fast1D, where the bc is respected
directly in integration).
"""
function dierckx_safe_integrate(
    spl::Dierckx.Spline1D,
    x1::FT,
    x2::FT;
    bc::BCT = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
) where {FT, BCT <: SSCF.Interpolation.ValidBoundaryConditions}

    y = FT(0)

    xp = FT.(extrema(spl.t)) # the outer knots define the end of the spline... this may or may not be the end of the data.
    yp = FT.(spl.(xp)) # the values at the knots

    spl_part = nothing
    bc_parts = nothing
    xbs = nothing
    fxbs = nothing
    dfdxbs = nothing
    if (x2 < xp[1])  # completely below
        spl_part = nothing
        bc_parts = ((x1, x2),)

        xbs = (xp[1],)
        fxbs = (yp[1],)
        if bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif bc isa SSCF.Interpolation.ErrorBoundaryCondition
            error("x2 < xp[1] and bc = $bc, this is not allowed")
        end

    elseif (x1 ≤ xp[1]) && (xp[1] ≤ x2 ≤ xp[end]) # partially below
        bc_parts = ((x1, xp[1]),)
        spl_part = (xp[1], x2)

        xbs = (xp[1],)
        fxbs = (yp[1],)
        if bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif (bc isa SSCF.Interpolation.ErrorBoundaryCondition) && (x1 < xp[1])
            error("x1 < xp[1] and bc = $bc, this is not allowed")
        end

    elseif (xp[1] ≤ x1) && (x2 ≤ xp[end]) # completely inside
        spl_part = (x1, x2)
        bc_parts = nothing
        return Dierckx.integrate(spl, x1, x2) # short circuit

    elseif (xp[1] ≤ x1 ≤ xp[end]) && (xp[end] ≤ x2) # partially above
        spl_part = (x1, xp[end])
        bc_parts = ((xp[end], x2),)

        xbs = (xp[end],)
        fxbs = (yp[end],)
        if (bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition)
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc isa SSCF.Interpolation.ErrorBoundaryCondition) && (x2 > xp[end])
            error("x2 > xp[end] and bc = $bc, this is not allowed")
        end

    elseif (xp[end] < x1) # completely above
        bc_parts = ((x1, x2),)
        spl_part = nothing

        xbs = (xp[end],)
        fxbs = (yp[end],)
        if (bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition)
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc isa SSCF.Interpolation.ErrorBoundaryCondition)
            error("x1 > xp[end] and bc = $bc, this is not allowed")
        end

    elseif (x1 < xp[1]) && (xp[end] < x2) # surrounding the data
        xbs = (xp[1], xp[end])
        fxbs = (yp[1], yp[end])

        spl_part = (xp[1], xp[end])
        bc_parts = ((x1, xp[1]), (xp[end], x2))

        if bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]), Dierckx.derivative(spl, xp[end]))
        elseif bc isa SSCF.Interpolation.ErrorBoundaryCondition
            error("x1 < xp[1] and x2 > xp[end] and bc = $bc, this is not allowed")
        end

    else
        error("Something went wrong with the interpolation")
    end

    if !isnothing(spl_part)
        y += Dierckx.integrate(spl, spl_part[1], spl_part[2])
    end

    if !isnothing(bc_parts)
        for (j, bc_part) in enumerate(bc_parts)
            a, b = bc_part
            if (bc isa SSCF.Interpolation.NearestBoundaryCondition)
                y += fxbs[j] * (b - a)
            elseif (bc isa SSCF.Interpolation.ExtrapolateBoundaryCondition)
                # linear expansion f(x) ≈ f(xb) + f'(xb)(x - xb) integrated over (a, b)
                y += (fxbs[j] * (b - a) + (dfdxbs[j] / 2) * ((b - xbs[j])^2 - (a - xbs[j])^2))
            end
        end
    end

    return y
end

SSCF.Interpolation.safe_integrate(spl::Dierckx.Spline1D, x1::FT, x2::FT; bc = SSCF.Interpolation.ExtrapolateBoundaryCondition()) where {FT} = dierckx_safe_integrate(spl, x1, x2; bc = bc)

end # module
