using ForwardDiff: ForwardDiff
using PCHIPInterpolation: PCHIPInterpolation
using Dierckx: Dierckx, Spline1D
# using Integrals: Integrals
using NonNegLeastSquares: NonNegLeastSquares, nonneg_lsq



# ============================================================================================================================================= #   

# Scalar fast linear interpolation
function fast1d_interp(
    x::TT,
    xp::AbstractVector{T},
    fp::AbstractVector{T};
    bc::BCT = ErrorBoundaryCondition(),
    monotonic::Bool = false,
    check_monotonic::Bool = false,
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}

    N = length(xp)

    # --- Optionally verify monotonicity (ascending OR descending) ---
    if monotonic
        asc = true
        desc = true
        if check_monotonic
            for i in 2:N
                if xp[i] < xp[i - 1]
                    asc = false
                end
                if xp[i] > xp[i - 1]
                    desc = false
                end
                if !asc && !desc
                    error("xp must be monotonic (ascending or descending)")
                end
            end
        else
            asc = xp[1] <= xp[end]
        end
    else
        asc = xp[1] <= xp[end]
    end

    # --- Out-of-bounds ---
    xmin, xmax = xp[1], xp[end]
    if asc
        if x < xmin
            if bc isa NearestBoundaryCondition
                return fp[1]
            elseif bc isa ExtrapolateBoundaryCondition
                return fp[1] + (fp[2] - fp[1]) * (x - xp[1]) / (xp[2] - xp[1])
            else
                error("x = $x below interpolation range [$xmin, $xmax]")
            end
        elseif x > xmax
            if bc isa NearestBoundaryCondition
                return fp[end]
            elseif bc isa ExtrapolateBoundaryCondition
                return fp[end - 1] + (fp[end] - fp[end - 1]) * (x - xp[end - 1]) / (xp[end] - xp[end - 1])
            else
                error("x = $x above interpolation range [$xmin, $xmax]")
            end
        end
    else
        # descending grid
        if x > xmin
            if bc isa NearestBoundaryCondition
                return fp[1]
            elseif bc isa ExtrapolateBoundaryCondition
                return fp[1] + (fp[2] - fp[1]) * (x - xp[1]) / (xp[2] - xp[1])
            else
                error("x = $x above interpolation range [$xmin, $xmax]")
            end
        elseif x < xmax
            if bc isa NearestBoundaryCondition
                return fp[end]
            elseif bc isa ExtrapolateBoundaryCondition
                return fp[end - 1] + (fp[end] - fp[end - 1]) * (x - xp[end - 1]) / (xp[end] - xp[end - 1])
            else
                error("x = $x below interpolation range [$xmin, $xmax]")
            end
        end
    end

    # --- Interval search ---
    if asc
        i = searchsortedlast(xp, x)
    else
        # descending: manual search without allocations
        i = 1
        while i < N && xp[i] > x
            i += 1
        end
        i -= 1  # interval is xp[i] >= x >= xp[i+1]
    end
    i = min(max(i, 1), N - 1)

    # --- Linear interpolation ---
    x0, x1 = xp[i], xp[i + 1]
    y0, y1 = fp[i], fp[i + 1]

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end





function fast1d_interp!(
    y::AbstractVector{TT},
    x::AbstractVector{TT},
    xp::AbstractVector{T},
    fp::AbstractVector{T};
    bc::BCT = ErrorBoundaryCondition(),
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}
    @inbounds for j in eachindex(x)
        y[j] = fast1d_interp(x[j], xp, fp; bc = bc)
    end
    return y
end

# Vector version (broadcastable)
function fast1d_interp(
    x::AbstractVector{TT},
    xp::AbstractVector{T},
    fp::AbstractVector{T};
    bc::BCT = ErrorBoundaryCondition(),
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}
    y = similar(x)
    fast1d_interp!(y, x, xp, fp; bc = bc)
    return y
end


struct Fast1DLinearInterpolant{X <: AbstractVector, BCT <: ValidBoundaryConditions}
    xp::X
    fp::X
    bc::BCT
end

"""
This constructor drops collinear points to optimize storage and performance
"""
function Fast1DLinearInterpolant(
    xp::X,
    fp::X;
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Bool = true,
) where {X <: AbstractVector, BCT <: ValidBoundaryConditions}
    @assert length(xp) == length(fp) "xp and fp must have the same length"

    if !drop_collinear || length(xp) <= 2
        return Fast1DLinearInterpolant{X, BCT}(xp, fp, bc)
    end

    # FT = eltype(X)
    # Always keep first point
    new_xp = FT[xp[1]]
    new_fp = FT[fp[1]]

    for i in 2:(length(xp) - 1)
        # Slopes
        slope1 = (fp[i] - fp[i - 1]) / (xp[i] - xp[i - 1])
        slope2 = (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
        if slope1 != slope2  # keep point if not collinear
            push!(new_xp, xp[i])
            push!(new_fp, fp[i])
        end
    end

    # Always keep last point
    push!(new_xp, xp[end])
    push!(new_fp, fp[end])

    # Convert back to original type
    new_xp = X(new_xp)
    new_fp = X(new_fp)
    return Fast1DLinearInterpolant{X, BCT}(new_xp, new_fp, bc)
end

""" Technically not type stable, so moved to separate method """
function Fast1DLinearInterpolant(
    xp::X,
    fp::X;
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Bool = true,
) where {X <: StaticArrays.StaticVector, BCT <: ValidBoundaryConditions}
    @assert length(xp) == length(fp) "xp and fp must have the same length"

    if !drop_collinear || length(xp) <= 2
        return Fast1DLinearInterpolant{X, BCT}(xp, fp, bc)
    end

    FT = eltype(X)

    # Always keep first point
    new_xp = FT[xp[1]]
    new_fp = FT[fp[1]]

    for i in 2:(length(xp) - 1)
        # Slopes
        slope1 = (fp[i] - fp[i - 1]) / (xp[i] - xp[i - 1])
        slope2 = (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
        if slope1 != slope2  # keep point if not collinear
            push!(new_xp, xp[i])
            push!(new_fp, fp[i])
        end
    end

    # Always keep last point
    push!(new_xp, xp[end])
    push!(new_fp, fp[end])

    # Convert back to original type
    out_type = StaticArrays.SVector{length(new_xp), eltype(X)}
    return Fast1DLinearInterpolant{out_type, BCT}(
        create_svector(new_xp)::out_type,
        create_svector(new_fp)::out_type,
        bc,
    )
end

# Make it broadcastable (for vector or static array evaluation)
Base.broadcastable(s::Fast1DLinearInterpolant) = Ref(s)

change_bc(
    sp::Fast1DLinearInterpolant{X, BCTO},
    new_bc::BCTN,
) where {X <: AbstractVector, BCTO <: ValidBoundaryConditions, BCTN <: ValidBoundaryConditions} =
    Fast1DLinearInterpolant{X, BCTN}(sp.xp, sp.fp, new_bc) # returns a new object, I guess you could use accessor or setfield or whatever
# Make it callable for scalar x
function (s::Fast1DLinearInterpolant)(x::T) where {T}
    fast1d_interp(x, s.xp, s.fp; bc = s.bc)
end

function fast1d_derivative(
    x::T,
    spl::Fast1DLinearInterpolant{X, BCT};
) where {T <: Real, X <: AbstractVector{T}, BCT <: ValidBoundaryConditions}
    xp, fp = spl.xp, spl.fp
    N = length(xp)

    bc = spl.bc

    if x < xp[1]
        if bc isa NearestBoundaryCondition
            return zero(T)
        elseif bc isa ExtrapolateBoundaryCondition
            return (fp[2] - fp[1]) / (xp[2] - xp[1])
        else
            error("x = $x below interpolation range and bc=$bc")
        end
    elseif x > xp[end]
        if bc isa NearestBoundaryCondition
            return zero(T)
        elseif bc isa ExtrapolateBoundaryCondition
            return (fp[end] - fp[end - 1]) / (xp[end] - xp[end - 1])
        else
            error("x = $x above interpolation range and bc=$bc")
        end
    else
        i = searchsortedlast(xp, x)
        if i == N
            i -= 1
        end
        return (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
    end
end


function fast1d_derivative!(
    y::AbstractVector{T}, # output array
    x::AbstractVector{T}, # derivative locations
    spl::Fast1DLinearInterpolant{X, BCT},
) where {T <: Real, X <: AbstractVector, BCT <: ValidBoundaryConditions}
    @inbounds for j in eachindex(x)
        y[j] = fast1d_derivative(x[j], spl)
    end
    return y
end
# Vectorized derivative
function fast1d_derivative(
    x::AbstractVector{T},
    spl::Fast1DLinearInterpolant{X, BCT},
) where {T <: Real, X <: AbstractVector, BCT <: ValidBoundaryConditions}
    y = similar(x)
    fast1d_derivative!(y, x, spl)
    return y
end



# ============================================================================================================================================= #   

# Build the spline for PCHIP
function build_spline(
    method::PCHIPInterpolationMethod,
    xp,
    fp;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    spl = PCHIPInterpolation.Interpolator(xp, fp)
    if bc isa ErrorBoundaryCondition
        error("Not implemented: pchip with bc=\"error\"")
    elseif bc isa ExtrapolateBoundaryCondition
        spl = pchip_extrapolate(spl)
    else
        error("Unsupported boundary condition type for PCHIP: $bc")
    end
    return spl
end

# Build the spline for PCHIPSmoothDerivative
function build_spline(
    method::PCHIPSmoothDerivativeInterpolationMethod,
    xp,
    fp;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    return pchip_smooth_derivative(
        xp,
        fp;
        bc = bc,
        f_enhancement_factor = method.f_enhancement_factor,
        f_p_enhancement_factor = method.f_p_enhancement_factor,
    )
end

# Build the spline for DierckxSpline1D
function build_spline(
    method::DierckxSpline1DInterpolationMethod,
    xp,
    fp;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    return Dierckx.Spline1D(xp, fp; k = method.k, bc = bc_string(bc))
end

# Build the spline for FastLinear1DInterpolationMethod
function build_spline(
    ::FastLinear1DInterpolationMethod,
    xp::X,
    fp::X;
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Bool = true,
) where {X <: AbstractVector, BCT <: ValidBoundaryConditions}
    # Return a closure that is broadcastable
    # return x -> fast1d_interp(x, xp, fp; bc=bc)
    return Fast1DLinearInterpolant(xp, fp; bc = bc, drop_collinear = drop_collinear)
end

# Build the spline for FastLinear1DInterpolationMethod
function build_spline(
    ::FastLinear1DInterpolationMethod,
    xp::AbstractVector{TX},
    fp::AbstractVector{TY};
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Bool = true,
) where {TX <: Real, TY <: Real, BCT <: ValidBoundaryConditions}
    # Return a closure that is broadcastable
    # return x -> fast1d_interp(x, xp, fp; bc=bc)
    FT = promote_type(TX, TY)
    return Fast1DLinearInterpolant(
        convert(Vector{FT}, xp),
        convert(Vector{FT}, fp);
        bc = bc,
        drop_collinear = drop_collinear,
    )
end

function build_spline(
    ::FastLinear1DInterpolationMethod,
    xp::StaticArrays.StaticVector{N, TX},
    fp::StaticArrays.StaticVector{N, TY};
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Bool = true,
) where {N, TX <: Real, TY <: Real, BCT <: ValidBoundaryConditions}
    # Return a closure that is broadcastable
    # return x -> fast1d_interp(x, xp, fp; bc=bc)
    FT = promote_type(TX, TY)
    return Fast1DLinearInterpolant(
        convert(StaticArrays.SVector{N, FT}, xp),
        convert(StaticArrays.SVector{N, FT}, fp);
        bc = bc,
        drop_collinear = drop_collinear,
    )
end

# Generic type-stable build_spline helper
function build_spline(
    ::Type{T},
    xp,
    fp;
    bc::BCT = ErrorBoundaryCondition(),
    k::Int = 1,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    drop_collinear::Bool = true,
) where {T <: AbstractInterpolationMethod, BCT <: ValidBoundaryConditions}
    if T <: PCHIPInterpolationMethod
        return build_spline(T(), xp, fp; bc = bc)
    elseif T <: PCHIPSmoothDerivativeInterpolationMethod
        concrete_method =
            T(f_enhancement_factor = f_enhancement_factor, f_p_enhancement_factor = f_p_enhancement_factor)
        return build_spline(concrete_method, xp, fp; bc = bc)
    elseif T <: DierckxSpline1DInterpolationMethod
        concrete_method = T(k)
        return build_spline(concrete_method, xp, fp; bc = bc)
    elseif T <: FastLinear1DInterpolationMethod
        return build_spline(T(), xp, fp; bc = bc, drop_collinear = drop_collinear)
    else
        error("Unknown interpolation method type: $T")
    end
end


# --------------------------------------------------------------------------------------------------------------------------------------------- #

# Generic interpolate function [[ Not type stable ]]
# Evaluate at x for each type
function interpolate_1d(
    x,
    xp,
    fp,
    method::PCHIPInterpolationMethod;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    spl = build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end

function interpolate_1d(
    x,
    xp,
    fp,
    method::PCHIPSmoothDerivativeInterpolationMethod;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    spl = build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end

function interpolate_1d(
    x,
    xp,
    fp,
    method::DierckxSpline1DInterpolationMethod;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    spl = build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end

function interpolate_1d(
    x,
    xp,
    fp,
    method::FastLinear1DInterpolationMethod;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    spl = build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end



function interpolate_1d(
    x,
    xp,
    fp;
    method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    bc::BCT = ErrorBoundaryCondition(),
    k::Int = 1,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    allow_Dierckx_k1_fastpath::Bool = true,
) where {BCT <: ValidBoundaryConditions}
    if method <: PCHIPInterpolationMethod
        return interpolate_1d(x, xp, fp, method(); bc = bc)
    elseif method <: PCHIPSmoothDerivativeInterpolationMethod
        concrete_method =
            method(f_enhancement_factor = f_enhancement_factor, f_p_enhancement_factor = f_p_enhancement_factor)
        return interpolate_1d(x, xp, fp, concrete_method; bc = bc)
    elseif method <: DierckxSpline1DInterpolationMethod
        if allow_Dierckx_k1_fastpath && (k == 1)
            # Fast linear path (should be fewer allocs than creating the Dierckx object
            return interpolate_1d(x, xp, fp, FastLinear1DInterpolation; bc = bc)
        end
        concrete_method = method(k)
        return interpolate_1d(x, xp, fp, concrete_method; bc = bc)
    elseif method <: FastLinear1DInterpolationMethod
        return interpolate_1d(x, xp, fp, method(); bc = bc)
    else
        error("Unknown interpolation method type: $method")
    end
end


# ============================================================================================================================================= #   

# TODO: Turn these piecewise pchip functions into a real function object (create my own class?) so we can analytically integrate them. it's just linear, quadratic and pchip

function pchip_extrapolate(spl::PCHIPInterpolation.Interpolator; return_piecewise_function::Bool = false)
    # extrapolating beyond the bounds of the data, we'll just continue the slope at the edge of the data for continuous derivative (though not necessarily smooth)
    xmin = spl.xs[1]
    xmax = spl.xs[end]
    ymin = spl.ys[1]
    ymax = spl.ys[end]

    return x -> begin
        if x < xmin
            dydx_xmin = PCHIPInterpolation._derivative(spl, Val(:begin), 1) # slope at beginning of interval 1
            return ymin - (xmin - x) * dydx_xmin
        elseif x > xmax
            dydx_xmax = PCHIPInterpolation._derivative(spl, Val(:end), length(spl.xs) - 1) # slope at end of interval N-1 (there are only N-1 intervals on N points) (N = length(spl.xs))
            return ymax + (x - xmax) * dydx_xmax
        else
            return spl(x)
        end
    end
end

"""
Create a function with a smooth second derivative by integrating a pchip approximation as the data's first derivative
Ideally you'd fit the spline to f(x) but f(x) may not be smooth. So we fit the spline to f'(x) and then integrate it to get f(x)
This comes at the cost of not necessarily hitting the points in f(x) exactly and rounding corners but it will be smooth...

To reduce the rounding, we can increase the resolution of the data by calculating the derivative spline at more points before creating the spline representation of the derivative
We do this via `enhancement_factor`, which allows us to add `enhancement_factor-1` points between each point in `xp` to constrain the spline of the spline derivative to more closely match the spline derivative and thus f(x) when integrated

We can also reduce rounding by boosting linear interpolation between the given data points, though this is somewhat of a guess as you don't actually know what the function looks like

Essentially:
    f_enhancement_factor: How closely the outcome matches linear interpolation between the data points. How rounded can the curve be?
    f_p_enhancement_factor: How closely the outcome matches the spline fit of f_p (possibly enhanced) -- i.e. does it cut corners or actually go to the points like the spline does?
# To Do: Support bc = NearestBoundaryCondition() and bc = FixedValueBoundaryCondition() for nearest neighbor and NaN respectively outside the bounds of xp
"""
function pchip_smooth_derivative(
    xp,
    fp;
    bc::BCT = ErrorBoundaryCondition(),
    # f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_enhancement_factor::Int = 1,
    # f_p_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Int = 1,
) where {BCT <: ValidBoundaryConditions}

    # if !isnothing(f_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
    if !isone(f_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
        # add enhancement_factor-1 points between each point (note this can lead to inexact errors if the interpolated points can't be cast to exterior type (like range needing float but x being in int))
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_enhancement_factor - 1))
        for i in 1:(length(xp) - 1)
            xp_new[((i - 1) * f_enhancement_factor + 1):(i * f_enhancement_factor)] .=
                range(xp[i], stop = xp[i + 1], length = f_enhancement_factor + 1)[1:(end - 1)]
        end
        xp_new[end] = xp[end]
        xp, fp = xp_new, interpolate_1d(xp_new, xp, fp; method = Fast1Linear1DInterpolationMethod, bc = bc)
    end

    # create a pchip interpolator to xp
    f_pchip_spl = PCHIPInterpolation.Interpolator(xp, fp) # should this be the extrapolatory one
    # differentiate it
    dfdx = ForwardDiff.derivative.(Ref(f_pchip_spl), xp)

    # if !isnothing(f_p_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
    if !isone(f_p_enhancement_factor)
        # add enhancement_factor-1 points between each point (note this can lead to inexact errors if the interpolated points can't be cast to exterior type (like range needing float but x being in int))
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_p_enhancement_factor - 1))
        for i in 1:(length(xp) - 1)
            xp_new[((i - 1) * f_p_enhancement_factor + 1):(i * f_p_enhancement_factor)] .=
                range(xp[i], stop = xp[i + 1], length = f_p_enhancement_factor + 1)[1:(end - 1)]
        end
        xp_new[end] = xp[end]
        xp, dfdx = xp_new, ForwardDiff.derivative.(Ref(f_pchip_spl), xp_new)
    end

    # create a pchip interpolator to that dxdp
    spl_dfdx = PCHIPInterpolation.Interpolator(xp, dfdx)
    dfp_dx_xmin = PCHIPInterpolation._derivative(spl_dfdx, Val(:begin), 1) # slope at beginning of interval 1
    dfp_dx_xmax = PCHIPInterpolation._derivative(spl_dfdx, Val(:end), length(xp) - 1) # slope at end of interval N-1 (there are only N-1 intervals on N points)
    dfdx_min = spl_dfdx.ys[1]
    dfdx_max = spl_dfdx.ys[end]
    ymin = fp[1]
    ymax = fp[end]

    xmin = xp[1]
    xmax = xp[end]

    # integrate it (only takes definite bounds so we can integrate between adjacent points and then cumsum?)
    return x -> begin
        if x < xp[1] # integrate from x to x_0
            if bc isa ErrorBoundaryCondition
                error(
                    "Requested x is below the minimum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate",
                )
            elseif bc isa ExtrapolateBoundaryCondition
                x_0 = xp[1] # aka xmin
                Δx = x_0 - x
                # assume derivative from x to x_0 is a line  f'(x) = x-> dfdx_min - dfp_dx_xmin * Δx, aka the second derivative is continous, the derivative is smooth | integrate to get f(x), but go from x_0 to x (so negative of x to x_0) 
                return ymin - (dfdx_min * Δx - dfp_dx_xmin * (+x_0 * Δx - (x_0^2 - x^2) / 2)) #+ (ymin+xmin)/dfdx_min + ymin  # -(...) bc we integrate from x_0 to x then take the negative 
            else
                error("Unsupported bc option $bc")
            end
        elseif x > xp[end] # integrate from x to x_N
            if bc isa ErrorBoundaryCondition
                error(
                    "Requested x is above the maximum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate",
                )
            elseif bc isa ExtrapolateBoundaryCondition
                x_0 = xp[end] # aka xmax
                Δx = x - x_0
                ymax = ymin + PCHIPInterpolation.integrate(spl_dfdx, xp[1], xp[end]) # more accurate bc integration means you're slightly off on the right side, not sure if it's just numerical error or what
                return ymax + (dfdx_max * Δx + dfp_dx_xmax * ((x^2 - x_0^2) / 2 - x_0 * Δx)) # -(...) bc we integrate from x_0 to x then take the negative 
            else
                error("Unsupported bc option $bc")
            end
        else
            return ymin + PCHIPInterpolation.integrate(spl_dfdx, xmin, x) # start at xmin,ymin and integrate to x
        end
    end

end




"""
1D Conservative Regridder

Given points (xp, yp), regrids to points (x,y) such that the integral of the function over the new points is (approximately) equal to the integral of the function over the old points
This is done by integrating the spline over the area of influence of each point, which is defined as the region over which the point is the nearest neighbor in the new grid. 

When upsampling, this should mostly just follow the spline (though slightly diffusive -- NOTE: is there any way to stop this diffusion from the "area of influence"? i.e. to be able to retrieve the original values? not really the spline is inherently diffusive of each point's delta fcn.
When downsampling though, you gain conservation when you otherwise would not have it.

TODO: Make a 2D version of this (including support for the 2D versio of conservative_spline_values for conservative regridding in 2D) (see `conservative_interpolate() in CalibrateEDMF.HelperFuncs`

"""
function conservative_regridder(
    x::AbstractVector,
    xp::AbstractVector,
    yp::AbstractVector;
    bc::BCT = ExtrapolateBoundaryCondition(), # must be extrapolate to integrate outside the knots (boundary of the data)... Dierckx integral ignored silently, we use dierckx_safe_integrate() instead
    k::Int = 1,
    method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    integrate_method::Symbol = :invert,
    rtol::FT = FT(1e-6),
    preserve_monotonicity::Bool = false,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::FT = FT(1e-8),
    enforce_conservation::Bool = true,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
) where {BCT <: ValidBoundaryConditions}

    # Build the spline
    # each node gets the mean of the spline over its area of influence, which we'll define as being the nearest neighbor
    # at the lower and upper end we'll extrapolate as far as the previous edge

    # get the spline
    spl = build_spline(
        method,
        xp,
        yp;
        bc = bc,
        k = k,
        f_enhancement_factor = f_enhancement_factor,
        f_p_enhancement_factor = f_p_enhancement_factor,
    )

    # calculate bin edges, at the end assume the same spacing, as is if it was a cell center
    n = length(x)
    x_edges = similar(x, n + 1)
    @inbounds begin
        x_edges[1] = x[1] - (x[2] - x[1]) / 2
        for i in 1:(n - 1)
            x_edges[i + 1] = (x[i] + x[i + 1]) / 2
        end
        x_edges[n + 1] = x[n] + (x[n] - x[n - 1]) / 2
    end



    y = zero.(x)

    for i in 1:length(x)
        # integrate over the area of influence

        if method <: DierckxSpline1DInterpolationMethod
            y[i] = dierckx_safe_integrate(spl, x_edges[i], x_edges[i + 1]; bc = bc) / (x_edges[i + 1] - x_edges[i])
        elseif method <: FastLinear1DInterpolationMethod
            y[i] = fast1d_safe_integrate(spl, x_edges[i], x_edges[i + 1]) / (x_edges[i + 1] - x_edges[i])
        elseif method isa AbstractPCHIPInterpolationMethod
            # we don't have a good way to interpolate the pchip output fcn bc it's piecewise, so we integrate it numerically
            # y[i] =
            #     Integrals.solve(
            #         Integrals.IntegralProblem((x, p) -> spl(x), (x_edges[i], x_edges[i + 1])),
            #         Integrals.QuadGKJL(),
            #     ).u / (x_edges[i + 1] - x_edges[i])
            # Remove dependence on Integrals for compat reasons with SciML (overly restrictive between Integrals 3 and 4). raise error
            error(
                "Conservative interpolation not yet supported for pchip due to compatibility issues with Integrals.jl vs 3.9, and 4.0, and other SciML packages like DiffEqBase.jl. Consider creating a version that relies on analytical solutions for extrapolation and pchip's integrate method inside.",
            )
        else
            error("method not recognized")
        end
    end


    if integrate_method == :integrate
        # Then we're done... 
    elseif integrate_method == :invert # This should be less diffusive and better preserve peaks so we make it the default
        # really, integrating the spline gives us the necessary mass, but it doesn't tell us the values at the points that would create a spline with that mass
        # really, what we want is a value such that calculating a new spline with the new points would yield the same total
        # for linear that's hard to do, but not so easy for a spline since you don't know what the new coefficients would be...

        # we need to find values that will yield the masses 

        # xf = x_edges
        # Δx = xf[2:end] .- xf[1:end-1]
        # spline_orig = spl
        # mc_avg = [first(Integrals.QuadGK.quadgk(spline_orig, xf[i], xf[i+1])) / Δx[i] for i in eachindex(Δx)]

        if enforce_positivity
            y_mean = mean(y)
            if (0 < y_mean < 1) # i think too large is fine, ignore small numbers, but bypass if mean is 0.
                # I think things start to break down around 2* eps(FT)^0.5, but we'll bring everything small up to 1 to be sure...
                y /= y_mean # we just scale up to 1, that should be fine, if the data has a huge range maybe something breaks but that's not a good fit for regridding like this anyway... [ maybe *= inv(y_mean) is faster?]
            end
        end

        _, y = conservative_spline_values(
            x_edges,
            y;
            k = k,
            method = method,
            bc = bc,
            f_enhancement_factor = f_enhancement_factor,
            f_p_enhancement_factor = f_p_enhancement_factor,
            rtol = rtol,
            enforce_positivity = enforce_positivity,
            nnls_alg = nnls_alg,
            nnls_tol = nnls_tol,
            A = A,
            Af = Af,
            yc = y,
        ) # this is the new y values that will yield the same mass as the original spline

        if enforce_positivity && (0 < y_mean < 1)
            y *= y_mean
        end

        if preserve_monotonicity # ensure our new values are within the bounds of their adjacent points in the original grid
            for i in 1:(length(x))


                # Find index of largest point in xp smaller than x_edges[i]
                i_low = searchsortedlast(xp, x_edges[i]) - (k - 1) # Find index of largest point in xp smaller than x_edges[i]. We've already a span of 1 from i_low to i_high so k-1 is enough width here I think... (meaning k=1, you'd still have at least 2 total points contributing... bc edges are all unique)
                # Find index of smallest point in xp larger than x_edges[i+1]
                i_high = searchsortedfirst(xp, x_edges[i + 1]) + (k - 1)  # Find index of smallest point in xp larger than x_edges[i+1]

                if (i_low < 1) || (i_high > length(xp)) # if our point has a support beyond the data bounds, we can't really hope to contain it... we could add an option to just only use the truncated data but...
                    continue
                end

                y[i] = clamp(y[i], minimum(yp[i_low:i_high]), maximum(yp[i_low:i_high])) # clamp the value to the bounds of the two points around it


            end
        end


        # inversion can be leaky, resolve if asked (:integrate is by definition not leaky since it's just the original integral of the spline.)
        if enforce_conservation
            if method <: DierckxSpline1DInterpolationMethod
                total = dierckx_safe_integrate(
                    Dierckx.Spline1D(x, y; k = k, bc = bc_string(bc)),
                    x_edges[1],
                    x_edges[end];
                    bc = bc,
                )
                if !iszero(total)
                    y *= dierckx_safe_integrate(spl, x_edges[1], x_edges[end]; bc = bc) / total # this is the total mass of the original data
                else
                    # we can't really do much here, maybe y is supposed to be nonzero but sum to zero, maybe not... who knows
                end
            elseif method <: FastLinear1DInterpolationMethod
                total = fast1d_safe_integrate(Fast1DLinearInterpolant(x, y; bc = bc), x_edges[1], x_edges[end])
                if !iszero(total)
                    y *= fast1d_safe_integrate(spl, x_edges[1], x_edges[end]) / total # this is the total mass of the original data
                else
                    # we can't really do much here, maybe y is supposed to be nonzero but sum to zero, maybe not... who knows
                end
            elseif method <: AbstractPCHIPInterpolationMethod
                # total =
                #     Integrals.solve(
                #         Integrals.IntegralProblem((x, p) -> spl(x), (x_edges[1], x_edges[end])),
                #         Integrals.QuadGKJL(),
                #     ).u
                error(
                    "Conservative interpolation not yet supported for pchip due to compatibility issues with Integrals.jl vs 3.9, and 4.0, and other SciML packages like DiffEqBase.jl. Consider creating a version that relies on analytical solutions for extrapolation and pchip's integrate method inside.",
                )
            end
        end


    else
        error("integrate_method not recognized")
        # each new point takes the average of its area of influence. we still use the spline in case we need to upsample
    end

    # this fcn doesn't really support returning a spline, you could take the new out and x and make a spline but it's kind of meaningless bc integrating that spline isn't conservative, it's the original spline that was the point...
    return y

end


# ----------------------------------------------------------------
# xc_from_xf! (in-place)
# ----------------------------------------------------------------
function xc_from_xf!(xc::AbstractVector{T}, xf::AbstractVector{T}) where {T <: AbstractFloat}
    N = length(xf)
    @assert length(xc) == N - 1
    @inbounds for i in 1:(N - 1)
        xc[i] = (xf[i] + xf[i + 1]) * T(0.5)
    end
    return xc
end

# ----------------------------------------------------------------
# xc_from_xf (allocating)
# ----------------------------------------------------------------
function xc_from_xf(xf::AbstractVector{T}) where {T <: AbstractFloat}
    N = length(xf)
    xc = similar(xf, N - 1)
    xc_from_xf!(xc, xf)
end

# ----------------------------------------------------------------
# xf_from_xc! (in-place)
# ----------------------------------------------------------------
function xf_from_xc!(xf::AbstractVector{T}, xc::AbstractVector{T}) where {T <: AbstractFloat}
    N = length(xc)
    @assert length(xf) == N + 1
    @inbounds begin
        for i in 1:(N - 1)
            xf[i + 1] = (xc[i] + xc[i + 1]) * T(0.5)
        end
        xf[1] = xc[1] - (xc[2] - xc[1]) * T(0.5)
        xf[N + 1] = xc[N] + (xc[N] - xc[N - 1]) * T(0.5)
    end
    return xf
end

# ----------------------------------------------------------------
# xf_from_xc (allocating)
# ----------------------------------------------------------------
function xf_from_xc(xc::AbstractVector{T}) where {T <: AbstractFloat}
    N = length(xc)
    xf = similar(xc, N + 1)
    xf_from_xc!(xf, xc)
end




"""
Given x_edges and integrated masses between them, calculates y at xc such that an equivalent spline on yc integrated over x_edges gives the same mass as the original spline.

E.g.

    x = [1,2,3,4,5,6]
    y = [0,0,1,-1,0,0]

    x_edges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

    spl = interpolate_1d(_, x, y; spl_kwargs... return_spl = true)
    masses = [first(Integrals.QuadGK.quadgk(spl, x_edges[i], x_edges[i+1])) for i in 1:(length(x_edges)-1)]

    x_new = [ 1.5, 3.5, 5.5]
    x_new_edges = [0.5, 2.5, 4.5, 6.5]
    We want to find y_new such integrating spl(_, x_new, y_new; spl_kwargs..., return_spl = true) over x_new_edges gives the same masses as the original spline...

    # Method
    Constructs a system `A * yc = mc`, where:
    - `m_total[i] = mc[i] * Δx[i]` is the total mass in cell `i`
    - `A[i,j] = ∫ φ_j(x) dx` over cell `i`, where `φ_j(x)` is a spline that is 1 at `xc[j]` and
    0 at all other centers `xc[k ≠ j]`.

    This approximates the effect of a spline basis function at each center, builds the mass
    conservation matrix, and solves the resulting linear system for the values `yc`.

    # Notes
    - This assumes `Dierckx.Spline1D` provides a usable approximation to a basis function by
    fitting a spline through a unit vector. Because Dierckx splines are global and not
    compactly supported, this method can be numerically unstable for large `k` or small `n`.
    - For robust conservative remapping, a local B-spline basis (e.g. from BSplineKit.jl)
    is recommended.


    Note this can go off the rails in extrapolation, and at high orders.
    It's also not that fast...
    

    Note -- `conservative_regridder calls this method...`

"""

function conservative_spline_values(
    xf::AbstractVector{FT},
    mc::AbstractVector{FT2},
    ;
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    return_spl::Bool = false,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    rtol::FT = FT(1e-6),
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::FT = FT(1e-8),
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
    yc::Union{AbstractVector{FT2}, Nothing} = nothing,
) where {FT <: Real, FT2 <: Real, BCT <: ValidBoundaryConditions}

    xc = xc_from_xf(xf)

    if (isnothing(Af) && !enforce_positivity) || (isnothing(A) && enforce_positivity) # we want to use Af because it's faster unless we are in enforce_positivity mode bc the NNLS solver can't use a factorization
        if isnothing(A)
            A = get_conservative_A(
                xc;
                bc = bc,
                k = k,
                method = method,
                f_enhancement_factor = f_enhancement_factor,
                f_p_enhancement_factor = f_p_enhancement_factor,
            )
        end
        if !enforce_positivity
            A = LinearAlgebra.factorize(A) # Replace A with its factorization for fast solves
        end
    else
        if !enforce_positivity # we know Af isn't nothing bc we checked above
            A = Af # just use the factorization if you have it
        end
    end

    if isnothing(yc)
        yc = zeros(FT2, length(xc)) # initialize yc if not given
    end

    if enforce_positivity
        if 0 < mean(mc) < (2 * eps(FT)^0.5)
            @warn "mean(mc) = $(mean(mc)) < [2 * eps($FT)^0.5 = $(2 * eps(FT)^0.5)]; this is very small, NNLS may arbitrarily converge to 0. Consider scaling your data to be larger."
        end
        yc .= max.(NonNegLeastSquares.nonneg_lsq(A, mc; alg = nnls_alg, tol = nnls_tol)[:], FT(0)) # Non-negative least square solver (bound the output by 0, since some algs like :pivot can leave underflow negatives like 1e-17 in NonNegativeLeastSquares.jl)
    else
        yc .= A \ mc # check 2nd run....
    end

    if any(!isfinite, mc)
        error("Received invalid (non-finite) input in mc = $mc")
    elseif any(!isfinite, yc) # && all(isfinite, mc) # moved condition to check above
        # @error("NaN values in yc from inputs: xf = $xf; mc = $mc; bc = $bc; A = $A; k = $k; return_spl = $return_spl; enforce_positivity = $enforce_positivity; nnls_alg = $nnls_alg")
        # set NaNs to zero  (is this good> we seemed to get from very very small numbers but idk..., like 1e-300
        resolve_nan!(yc)
    end
    return xc, yc
end

function contiguous_true_ranges(bits::BitVector)
    len = length(bits)
    ranges = Vector{UnitRange{Int}}()

    i = 1
    @inbounds while i <= len
        if bits[i]
            start = i
            i += 1
            while i <= len && bits[i]
                i += 1
            end
            push!(ranges, start:(i - 1))
        else
            i += 1
        end
    end

    return ranges
end

function contiguous_true_ranges(mat::BitArray{2}; dim::Int = 1)
    # build a 1D Bool mask of length = size(mat, dim)
    mask = if dim == 1
        BitVector(any(mat[i, :]) for i in 1:size(mat, 1)) # in each row, any true [ so output is along the column ]
    elseif dim == 2
        BitVector(any(mat[:, j]) for j in 1:size(mat, 2)) # in each column, any true [ so output is along the row ]
    else
        throw(ArgumentError("`dim` must be 1 (rows) or 2 (cols), got $dim"))
    end

    # delegate to the vector version
    contiguous_true_ranges(mask)
end

"""
    2D version
"""
function conservative_spline_values(
    xf::AbstractVector{FT},
    mc::AbstractMatrix{FT}, # this is the mass concentration, either a vector or a matrix
    ;
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    return_spl::Bool = false,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    rtol::FT = FT(1e-6),
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::FT = FT(1e-8), # default for Float64 in package
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
    yc::Union{AbstractMatrix{FT2}, Nothing} = nothing,
    inplace::Bool = false, # if true, we modify mc inplace and return nothing, otherwise we return the modified mc
) where {FT <: Real, FT2 <: Real, BCT <: ValidBoundaryConditions}

    n = length(xf) - 1
    nx = n
    xc = xc_from_xf(xf)

    squeeze_y::Bool = false # if mc is a vector, we will squeeze it back to a vector later
    if ndims(mc) == 1
        mc = reshape(mc, n, 1) # make it a column vector if it's a vector [does this mess w/ inplace? idk...] doesn't reassign any memory, iguess you can pass in mc as a view if you really want to be inplace
        squeeze_y = true # we'll squeeze it back later
    end
    nt = size(mc, 2) # number of time steps, we assume that the second dimension is the time dimension


    # Get A
    if (isnothing(Af) && !enforce_positivity) || (isnothing(A) && enforce_positivity) # we want to use Af because it's faster unless we are in enforce_positivity mode bc the NNLS solver can't use a factorization
        if isnothing(A)
            A = get_conservative_A(
                xc;
                bc = bc,
                k = k,
                method = method,
                f_enhancement_factor = f_enhancement_factor,
                f_p_enhancement_factor = f_p_enhancement_factor,
            )
        end
        if !enforce_positivity
            A = LinearAlgebra.factorize(A) # Replace A with its factorization for fast solves
        end
    else
        if !enforce_positivity # we know Af isn't nothing bc we checked above
            A = Af # just use the factorization if you have it
        end
    end


    # Do interpolation
    if inplace
        yc = mc # modify mc inplace
    else
        yc = zeros(FT, nx, nt) # initialize yc as a zero matrix
    end

    # @time if enforce_positivity
    if enforce_positivity # if we want to enforce positivity, we need to use a non-negative least squares solver
        # if 0 < mean(mc) < (2 * eps(FT)^0.5)
        # @info "mc = $(mc); mean(mc) = $(mean(mc)); eps(FT) = $(eps(FT)); 2 * eps(FT)^0.5 = $(2 * eps(FT)^0.5)"
        if any(x -> (0 < x < 2 * eps(FT)^0.5), mean(mc, dims = 1)) # check if the mean of mc is small, if so, we should warn the user
            @warn "mean(mc) = $(mean(mc)) < [2 * eps($FT)^0.5 = $(2 * eps(FT)^0.5)]; this is very small, NNLS may arbitarily converge to 0. Consider scaling your data to be larger."
        end

        # whereever M is 0, we dont need to solve there
        # Probably faster for things like profiles of condensate w/ many all zero rows. qt and theta it wont help though (but we don't really use those in outputs, only in postprocessing.)

        valid_inds = (!iszero).(mc) # this is a boolean vector of where mc is nonzero

        # we could still do all the solves in 1 and it seemed to save about 33% but idk how that plays w/ contiguous regions...

        for contiguous_region in contiguous_true_ranges(valid_inds; dim = 1) # any valid in time we keep, so check each row [dim = 1 is left behind as one column]
            A_sub = @view A[contiguous_region, contiguous_region] # get the submatrix of A for the contiguous region
            mc_sub = @view mc[contiguous_region, :] # get the subvector of mc for the contiguous region

            yc[contiguous_region, :] .=
                max.(NonNegLeastSquares.nonneg_lsq(A_sub, mc_sub; alg = nnls_alg, tol = nnls_tol), FT(0)) # Non-negative least square solver  (bound the output by 0, since some algs like :pivot can leave underflow negatives like 1e-17 in NonNegativeLeastSquares.jl). As of yet there's no inplace NNLS so we keep the braodcast w/ no preallocation
        end

    else
        yc .= A \ mc # solve for yc using the factorization of A
    end

    if any(!isfinite, mc)
        error(
            "Received invalid (non-finite) input in mc = $mc. NaNs inputs are not supported for conservative regridding because they break the matrix solve.",
        )
    elseif any(!isfinite, yc) # && all(isfinite, mc) # moved condition to check above
        @error(
            "NaN values in yc from inputs: xf = $xf; mc = $mc; bc = $bc; A = $A; k = $k; return_spl = $return_spl; enforce_positivity = $enforce_positivity; nnls_alg = $nnls_alg"
        )
        # set NaNs to zero  (is this good> we seemed to get from very very small numbers but idk..., like 1e-300
        resolve_nan!(yc)
    end

    if squeeze_y
        yc = vec(yc) # squeeze the output back to a vector if it was a vector (i.e. undo any reshape)
    end

    return xc, yc

end

function get_conservative_A(
    xc::AbstractVector{FT};
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
) where {FT, BCT <: ValidBoundaryConditions}

    # calculate bin edges, at the end assume the same spacing, as is if it was a cell center
    n = length(xc)
    xf = similar(xc, n + 1)
    @inbounds begin
        xf[1] = xc[1] - (xc[2] - xc[1]) / 2
        for i in 1:(n - 1)
            xf[i + 1] = (xc[i] + xc[i + 1]) / 2
        end
        xf[n + 1] = xc[n] + (xc[n] - xc[n - 1]) / 2
    end

    n = length(xc)
    if FT <: Int
        A = zeros(Float64, n, n) # probably cant factorize/integrate/transform in int... need a float.
    else
        A = zeros(FT, n, n)
    end

    for j in 1:n
        ej = zeros(FT, n)
        ej[j] = FT(1.0)
        φj = build_spline(
            method,
            xc,
            ej;
            k = k,
            bc = bc,
            f_enhancement_factor = f_enhancement_factor,
            f_p_enhancement_factor = f_p_enhancement_factor,
        )


        if method ∈ (DierckxSpline1DInterpolationMethod, FastLinear1DInterpolationMethod)
            width = k + 1 # Support width of each spline basis fcn [ this should push us from quadratic n^2 to linear kn ]
            i_start = max(1, j - width)
            i_end = min(n, j + width)
            for i in i_start:i_end
                # A[i, j] = safe_integrate(φj, xf[i], xf[i + 1]; bc = bc) / (xf[i + 1] - xf[i])
                if method <: FastLinear1DInterpolationMethod
                    A[i, j] = fast1d_safe_integrate(φj, xf[i], xf[i + 1]) / (xf[i + 1] - xf[i]) # this is fast
                else
                    A[i, j] = dierckx_safe_integrate(φj, xf[i], xf[i + 1]; bc = bc) / (xf[i + 1] - xf[i]) # this is fast
                end
            end

        elseif method <: AbstractPCHIPInterpolationMethod
            # for i in 1:n
            #     A[i, j] = first(Integrals.QuadGK.quadgk(φj, xf[i], xf[i + 1])) / (xf[i + 1] - xf[i]) # it could be anything so integrate (note this is very slow...)
            # end
            error(
                "Conservative interpolation not yet supported for pchip due to compatibility issues with Integrals.jl vs 3.9, and 4.0, and other SciML packages like DiffEqBase.jl. Consider creating a version that relies on analytical solutions for extrapolation and pchip's integrate method inside.",
            )
        end
    end
    A .= max.(A, FT(0)) # ensure no negative values, sometimes we get like tiny 1e-17s, not sure why, might be like subtraction underflow
    return A
end







"""
    rn Dierckx doesn't respect boundary conditions mode (extrapolate, nearest etc) when integrating, so we have to do it ourselves.
    This differs from fast1d where the bc is respected in integration.
"""
function dierckx_safe_integrate(
    spl::Dierckx.Spline1D,
    x1::FT,
    x2::FT,
    ;
    bc::BCT = ExtrapolateBoundaryCondition(),
) where {FT, BCT <: ValidBoundaryConditions}

    y = FT(0)

    xp = FT.(extrema(spl.t)) # the outer knots define the end of the spline... this may or may not be the end of the data.
    yp = FT.(spl.(xp)) # the values at the knots


    spl_part = nothing
    bc_parts = nothing
    xbs = nothing
    fxbs = nothing
    dfdxbs = Nothing
    if (x2 < xp[1])  # were completely below
        spl_part = nothing
        bc_parts = ((x1, x2),)

        xbs = (xp[1],)
        fxbs = (yp[1],)
        if bc isa ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif bc isa ErrorBoundaryCondition
            error("x2 < xp[1] and bc = $bc, this is not allowed")
        end

    elseif (x1 ≤ xp[1]) && (xp[1] ≤ x2 ≤ xp[end]) # we're only partially outside the bounds of the data
        bc_parts = ((x1, xp[1]),)
        spl_part = (xp[1], x2)

        xbs = (xp[1],)
        fxbs = (yp[1],)
        if bc isa ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif (bc isa ErrorBoundaryCondition) && (x1 < xp[1])
            error("x1 < xp[1] and bc = $bc, this is not allowed")
        end

    elseif (xp[1] ≤ x1) && (x2 ≤ xp[end]) # were completely inside the bounds of the data
        spl_part = (x1, x2)
        bc_parts = nothing
        return Dierckx.integrate(spl, x1, x2) # short circuit

    elseif (xp[1] ≤ x1 ≤ xp[end]) && (xp[end] ≤ x2) # we're only partially outside the bounds of the data
        spl_part = (x1, xp[end])
        bc_parts = ((xp[end], x2),)

        xbs = (xp[end],)
        fxbs = (yp[end],)
        if (bc isa ExtrapolateBoundaryCondition)
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc isa ErrorBoundaryCondition) && (x2 > xp[end])
            error("x2 > xp[end] and bc = $bc, this is not allowed")
        end

    elseif (xp[end] < x1) # were completely above the bounds of the data
        bc_parts = ((x1, x2),)
        spl_part = nothing

        xbs = (xp[end],)
        fxbs = (yp[end],)
        if (bc isa ExtrapolateBoundaryCondition)
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc isa ErrorBoundaryCondition)
            error("x1 > xp[end] and bc = $bc, this is not allowed")
        end


    elseif (x1 < xp[1]) && (xp[end] < x2) # we completely surround the data

        xbs = (xp[1], xp[end])
        fxbs = (yp[1], yp[end])

        spl_part = (xp[1], xp[end])
        bc_parts = ((x1, xp[1]), (xp[end], x2))

        if bc isa ExtrapolateBoundaryCondition
            dfdxbs = (Dierckx.derivative(spl, xp[1]), Dierckx.derivative(spl, xp[end]))
        elseif bc isa ErrorBoundaryCondition
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
            if (bc isa NearestBoundaryCondition)
                # This is just a constant "

                y += fxbs[j] * (b - a)
            elseif (bc isa ExtrapolateBoundaryCondition)
                # calculate the integral using derivative and values at the edge
                #=
                    Approximate f(x) using a linear expansion around a known point xb:
                        f(x) ≈ fxb + dfdxb * (x - xb)
                        where:
                            f(xb) = fxb
                            f'(xb) = dfdxb

                    Then the antiderivative y(x) = ∫ f(x) dx is:
                        y(x) = fxb * (x - xb) + (dfdxb / 2) * (x - xb)^2

                    To compute the integral of f(x) over the interval (a, b) = (bc_part[1], bc_part[2]):
                        ∫[a to b] f(x) dx ≈ [fxb * (x - xb) + (dfdxb / 2) * (x - xb)^2] evaluated from x = a to x = b

                        So:
                        ∫[a to b] f(x) dx ≈ [fxb*(b - xb) + (dfdxb/2)*(b - xb)^2] - [fxb*(a - xb) + (dfdxb/2)*(a - xb)^2]
                                        = fxb * (b - a) + (dfdxb / 2) * [ (b - xb)^2 - (a - xb)^2 ]
                        The average is just all that divided by (b - a)
                =#
                y += (fxbs[j] * (b - a) + (dfdxbs[j] / 2) * ((b - xbs[j])^2 - (a - xbs[j])^2))


            end
        end
    end

    return y
end


"""
Safe integration of Fast1DLinearInterpolant, respecting bc modes.

Unlike Dierckx, Fast1D respects boundary conditions in integration, so we do not pass bc in
"""
function fast1d_safe_integrate(
    spl::Fast1DLinearInterpolant{X, BCT},
    x1::T1,
    x2::T2;
) where {T1 <: Real, T2 <: Real, X <: AbstractVector, BCT <: ValidBoundaryConditions}
    ST = eltype(spl.xp)
    FT = promote_type(T1, T2, ST)

    bc = spl.bc

    y = zero(FT)
    xp, fp = spl.xp, spl.fp
    N = length(xp)

    # Internal derivative function
    deriv(x) = fast1d_derivative(x, spl)

    # # Prepare BCs as tuples
    spl_part = nothing

    bc_parts::Tuple{Vararg{Tuple{FT, FT}}} = ()
    xbs::Tuple{Vararg{FT}} = ()
    fxbs::Tuple{Vararg{FT}} = ()
    dfdxbs::Tuple{Vararg{FT}} = ()


    if x2 < xp[1]  # completely below
        bc_parts = ((x1, x2),)
        xbs = (xp[1],)
        fxbs = (fp[1],)
        dfdxbs = (bc isa ExtrapolateBoundaryCondition) ? (deriv(xp[1]),) : dfdxbs
        (bc isa ErrorBoundaryCondition) && error("x2 < xp[1] and bc=$bc not allowed")

    elseif x1 <= xp[1] <= x2 <= xp[end]  # partially below
        spl_part = (xp[1], x2)
        bc_parts = ((x1, xp[1]),)
        xbs = (xp[1],)
        fxbs = (fp[1],)
        dfdxbs = (bc isa ExtrapolateBoundaryCondition) ? (deriv(xp[1]),) : dfdxbs
        (bc isa ErrorBoundaryCondition) && (x1 < xp[1]) && error("x1 < xp[1] and bc=$bc not allowed")

    elseif xp[1] <= x1 && x2 <= xp[end]  # completely inside
        return fast1d_integrate_internal(spl, x1, x2)

    elseif xp[1] <= x1 <= xp[end] <= x2  # partially above
        spl_part = (x1, xp[end])
        bc_parts = ((xp[end], x2),)
        xbs = (xp[end],)
        fxbs = (fp[end],)
        dfdxbs = (bc isa ExtrapolateBoundaryCondition) ? (deriv(xp[end]),) : dfdxbs
        (bc isa ErrorBoundaryCondition) && (x2 > xp[end]) && error("x2 > xp[end] and bc=$bc not allowed")

    elseif xp[end] < x1  # completely above
        bc_parts = ((x1, x2),)
        xbs = (xp[end],)
        fxbs = (fp[end],)
        dfdxbs = (bc isa ExtrapolateBoundaryCondition) ? (deriv(xp[end]),) : dfdxbs
        (bc isa ErrorBoundaryCondition) && error("x1 > xp[end] and bc=$bc not allowed")

    elseif x1 < xp[1] && xp[end] < x2  # surrounding
        spl_part = (xp[1], xp[end])
        bc_parts = ((x1, xp[1]), (xp[end], x2))
        xbs = (xp[1], xp[end])
        fxbs = (fp[1], fp[end])
        dfdxbs = (bc isa ExtrapolateBoundaryCondition) ? (deriv(xp[1]), deriv(xp[end])) : dfdxbs
        (bc isa ErrorBoundaryCondition) && error("x1 < xp[1] and x2 > xp[end] and bc=$bc not allowed")
    end

    # integrate interior part
    if !isnothing(spl_part)
        y += fast1d_integrate_internal(spl, spl_part[1], spl_part[2])
    end

    # integrate bc parts
    # if !isnothing(bc_parts)
    if !isempty(bc_parts)
        for j in eachindex(bc_parts)
            a, b = bc_parts[j]
            if bc isa NearestBoundaryCondition
                y += fxbs[j] * (b - a)
            elseif bc isa ExtrapolateBoundaryCondition
                y += fxbs[j] * (b - a) + (dfdxbs[j] / 2) * ((b - xbs[j])^2 - (a - xbs[j])^2)
            end
        end
    end

    return y
end


# Internal integration assuming fully inside xp range
function fast1d_integrate_internal(
    spl::Fast1DLinearInterpolant{X, BCT},
    x1::T1,
    x2::T2,
) where {T1 <: Real, T2 <: Real, X <: AbstractVector, BCT <: ValidBoundaryConditions}
    ST = eltype(spl.xp)
    FT = promote_type(T1, T2, ST)

    xp, fp = spl.xp, spl.fp
    N = length(xp)
    y = zero(FT)

    # Loop over each interval intersecting [x1,x2]
    i1 = searchsortedlast(xp, x1)
    i2 = searchsortedlast(xp, x2)
    i1 = clamp(i1, 1, N - 1)
    i2 = clamp(i2, 1, N - 1)

    for i in i1:i2
        xa = max(x1, xp[i])
        xb = min(x2, xp[i + 1])
        slope = (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
        y += fp[i] * (xb - xa) + slope / 2 * (xb - xa)^2
    end

    return y
end

"""
BC is only passed for Dierckx splines, Fast1D respects bc in integration already.
"""
safe_integrate(spl, x1::FT, x2::FT; bc = ExtrapolateBoundaryCondition()) where {FT} =
    error("safe_integrate not implemented for spline type $(typeof(spl))")
safe_integrate(spl::Dierckx.Spline1D, x1::FT, x2::FT; bc = ExtrapolateBoundaryCondition()) where {FT} =
    dierckx_safe_integrate(spl, x1, x2; bc = bc)
safe_integrate(
    spl::Fast1DLinearInterpolant{X, BCT},
    x1::FT,
    x2::FT;
    bc = ExtrapolateBoundaryCondition(),
) where {X <: AbstractVector, BCT <: ValidBoundaryConditions, FT} = fast1d_safe_integrate(spl, x1, x2) # We do not pass bc through
