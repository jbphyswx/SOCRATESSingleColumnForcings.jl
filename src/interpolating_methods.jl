using ForwardDiff: ForwardDiff
using PCHIPInterpolation: PCHIPInterpolation
using Dierckx: Dierckx, Spline1D
# using Integrals: Integrals
using NonNegLeastSquares: NonNegLeastSquares, nonneg_lsq




function pyinterp(
    x,
    xp,
    fp;
    bc::String = "error",
    k::Int = 1,
    method::Symbol = :Spline1D,
    return_spl::Bool = false,
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
)
    if method == :pchip
        spl = PCHIPInterpolation.Interpolator(xp, fp) # Has a continuous (but not smooth) derivative (piecewise cubic hermite interpolating polynomial, is unique so no k)
        if bc == "error"

        elseif bc == "extrapolate"
            spl = pchip_extrapolate(spl) # allow extrapolating beyond the edges
        end

    elseif method ∈ [:pchip_smooth_derivative, :pchip_smooth]
        spl = pchip_smooth_derivative(
            xp,
            fp;
            bc = bc,
            f_enhancement_factor = f_enhancement_factor,
            f_p_enhancement_factor = f_p_enhancement_factor,
        ) # how our current implementation works
    elseif method ∈ [:Spline1D, :Dierckx]
        spl = Dierckx.Spline1D(xp, fp; k = k, bc = bc) # k=1 has disjoint derivative, k=2,3 is smooth in derivative but not monotonicity preserving... this is matlab pchip
    else
        error("method not recognized")
    end
    if return_spl
        return spl
    else
        return spl.(vec(x)) # have to broadcast bc pchip requires it...
    end
end

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
# To Do: Support bc = "nearest" and bc = "NaN" for nearest neighbor and NaN respectively outside the bounds of xp
"""
function pchip_smooth_derivative(
    xp,
    fp;
    bc::String = "error",
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
)

    if !isnothing(f_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
        # add enhancement_factor-1 points between each point (note this can lead to inexact errors if the interpolated points can't be cast to exterior type (like range needing float but x being in int))
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_enhancement_factor - 1))
        for i in 1:(length(xp) - 1)
            xp_new[((i - 1) * f_enhancement_factor + 1):(i * f_enhancement_factor)] .=
                range(xp[i], stop = xp[i + 1], length = f_enhancement_factor + 1)[1:(end - 1)]
        end
        xp_new[end] = xp[end]
        xp, fp = xp_new, pyinterp(xp_new, xp, fp; method = :Dierckx, return_spl = false)
    end

    # create a pchip interpolator to xp
    f_pchip_spl = PCHIPInterpolation.Interpolator(xp, fp) # should this be the extrapolatory one
    # differentiate it
    dfdx = ForwardDiff.derivative.(Ref(f_pchip_spl), xp)

    if !isnothing(f_p_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
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
            if bc == "error"
                error(
                    "Requested x is below the minimum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate",
                )
            elseif bc == "extrapolate"
                x_0 = xp[1] # aka xmin
                Δx = x_0 - x
                # assume derivative from x to x_0 is a line  f'(x) = x-> dfdx_min - dfp_dx_xmin * Δx, aka the second derivative is continous, the derivative is smooth | integrate to get f(x), but go from x_0 to x (so negative of x to x_0) 
                return ymin - (dfdx_min * Δx - dfp_dx_xmin * (+x_0 * Δx - (x_0^2 - x^2) / 2)) #+ (ymin+xmin)/dfdx_min + ymin  # -(...) bc we integrate from x_0 to x then take the negative 
            else
                error("Unsupported bc option $bc")
            end
        elseif x > xp[end] # integrate from x to x_N
            if bc == "error"
                error(
                    "Requested x is above the maximum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate",
                )
            elseif bc == "extrapolate"
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

"""
function conservative_regridder(
    x::AbstractVector,
    xp::AbstractVector,
    yp::AbstractVector;
    bc::String = "extrapolate", # must be extrapolate to integrate outside the knots (boundary of the data)... Dierckx integral ignored silently, we use dierckx_safe_integrate() instead
    k::Int = 1,
    method::Symbol = :Spline1D,
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
    integrate_method::Symbol = :invert,
    rtol = FT(1e-6),
    preserve_monotonicity::Bool = false,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol = FT(1e-8),
    enforce_conservation::Bool = true,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, Nothing} = nothing,
) where {}

    # Build the spline
    # each node gets the mean of the spline over its area of influence, which we'll define as being the nearest neighbor
    # at the lower and upper end we'll extrapolate as far as the previous edge

    # get the spline
    spl = pyinterp(
        x,
        xp,
        yp;
        k = k,
        bc = bc,
        method = method,
        f_enhancement_factor = f_enhancement_factor,
        f_p_enhancement_factor = f_p_enhancement_factor,
        return_spl = true,
    )

    # calculate bin edges, at the end assume the same spacing, as is if it was a cell center
    n = length(x)
    x_edges = Vector{FT}(undef, length(x) + 1)
    @inbounds begin
        x_edges[1] = x[1] - (x[2] - x[1]) / 2
        for i in 1:(n - 1)
            x_edges[i + 1] = (x[i] + x[i + 1]) / 2
        end
        x_edges[n + 1] = x[n] + (x[n] - x[n - 1]) / 2
    end



    y::Vector{FT} = zeros(FT, size(x))
    for i in 1:length(x)
        # integrate over the area of influence

        if method == :Spline1D
            y[i] = dierckx_safe_integrate(spl, x_edges[i], x_edges[i + 1]; bc = bc) / (x_edges[i + 1] - x_edges[i])
        elseif method ∈ (:pchip, :pchip_smooth_derivative, :pchip_smooth)
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
            if method == :Spline1D
                total =
                    dierckx_safe_integrate(Dierckx.Spline1D(x, y; k = k, bc = bc), x_edges[1], x_edges[end]; bc = bc)
                if !iszero(total)
                    y *= dierckx_safe_integrate(spl, x_edges[1], x_edges[end]; bc = bc) / total # this is the total mass of the original data
                else
                    # we can't really do much here, maybe y is supposed to be nonzero but sum to zero, maybe not... who knows
                end
            elseif method ∈ (:pchip, :pchip_smooth_derivative, :pchip_smooth)
                # total =
                #     Integrals.solve(
                #         Integrals.IntegralProblem((x, p) -> spl(x), (x_edges[1], x_edges[end])),
                #         Integrals.QuadGKJL(),
                #     ).u
                error(
                    "Conservative interpolation not yet supported for pchip due to compatibility issues with Integrals.jl vs 3.9, and 4.0, and other SciML packages like DiffEqBase.jl. Consider creating a version that relies on analytical solutions for extrapolation and pchip's integrate method inside.",
                )
                if !iszero(total)
                    y *= dierckx_safe_integrate(spl, x_edges[1], x_edges[end]; bc = bc) / total # this is the total mass of the original data
                else
                    # we can't really do much here, maybe y is supposed to be nonzero but sum to zero, maybe not... who knows
                end
            end
        end


    else
        error("integrate_method not recognized")
        # each new point takes the average of its area of influence. we still use the spline in case we need to upsample
    end

    # this fcn doesn't really support returning a spline, you could take the new out and x and make a spline but it's kind of meaningless bc integrating that spline isn't conservative, it's the original spline that was the point...
    return y

end


"""
Given x_edges and integrated masses between them, calculates y at xc such that an equivalent spline on yc integrated over x_edges gives the same mass as the original spline.

E.g.

    x = [1,2,3,4,5,6]
    y = [0,0,1,-1,0,0]

    x_edges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

    spl = pyinterp(_, x, y; spl_kwargs... return_spl = true)
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
    

"""

function conservative_spline_values(
    xf::AbstractVector{FT},
    mc::AbstractVector{FT2},
    ;
    bc::String = "extrapolate",
    k::Int = 1,
    method::Symbol = :Spline1D,
    return_spl::Bool = false,
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
    rtol = FT(1e-6),
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol = FT(1e-8),
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, Nothing} = nothing,
    yc::Union{AbstractVector{FT2}, Nothing} = nothing,
) where {FT, FT2}

    xc = FT(0.5) .* (xf[1:(end - 1)] .+ xf[2:end])

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


function get_conservative_A(
    xc::AbstractVector{FT};
    bc::String = "extrapolate",
    k::Int = 1,
    method::Symbol = :Spline1D,
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
) where {FT}

    # calculate bin edges, at the end assume the same spacing, as is if it was a cell center
    n = length(xc)
    xf = Vector{FT}(undef, length(xc) + 1)
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
        φj = pyinterp(
            xc, # we're not using the avlue
            xc,
            ej;
            k = k,
            bc = bc,
            method = :Spline1D,
            return_spl = true, # we want a spline out
            f_enhancement_factor = f_enhancement_factor,
            f_p_enhancement_factor = f_p_enhancement_factor,
        )


        if method == :Spline1D
            width = k + 1 # Support width of each spline basis fcn [ this should push us from quadratic n^2 to linear kn ]
            i_start = max(1, j - width)
            i_end = min(n, j + width)
            for i in i_start:i_end
                A[i, j] = dierckx_safe_integrate(φj, xf[i], xf[i + 1]; bc = bc) / (xf[i + 1] - xf[i]) # this is fast
            end

        elseif method ∈ (:pchip, :pchip_smooth_derivative, :pchip_smooth)
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
"""
function dierckx_safe_integrate(spl::Dierckx.Spline1D, x1::FT, x2::FT, ; bc::String = "extrapolate") where {FT}

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
        if bc == "extrapolate"
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif bc == "error"
            error("x2 < xp[1] and bc = $bc, this is not allowed")
        end

    elseif (x1 ≤ xp[1]) && (xp[1] ≤ x2 ≤ xp[end]) # we're only partially outside the bounds of the data
        bc_parts = ((x1, xp[1]),)
        spl_part = (xp[1], x2)

        xbs = (xp[1],)
        fxbs = (yp[1],)
        if bc == "extrapolate"
            dfdxbs = (Dierckx.derivative(spl, xp[1]),)
        elseif (bc == "error") && (x1 < xp[1])
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
        if bc == "extrapolate"
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc == "error") && (x2 > xp[end])
            error("x2 > xp[end] and bc = $bc, this is not allowed")
        end

    elseif (xp[end] < x1) # were completely above the bounds of the data
        bc_parts = ((x1, x2),)
        spl_part = nothing

        xbs = (xp[end],)
        fxbs = (yp[end],)
        if bc == "extrapolate"
            dfdxbs = (Dierckx.derivative(spl, xp[end]),)
        elseif (bc == "error")
            error("x1 > xp[end] and bc = $bc, this is not allowed")
        end


    elseif (x1 < xp[1]) && (xp[end] < x2) # we completely surround the data

        xbs = (xp[1], xp[end])
        fxbs = (yp[1], yp[end])

        spl_part = (xp[1], xp[end])
        bc_parts = ((x1, xp[1]), (xp[end], x2))

        if bc == "extrapolate"
            dfdxbs = (Dierckx.derivative(spl, xp[1]), Dierckx.derivative(spl, xp[end]))
        elseif bc == "error"
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
            if bc == "nearest"
                # This is just a constant "

                y += fxbs[j] * (b - a)
            elseif bc == "extrapolate"
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

                # elseif bc == "zero"
                #     # do nothing
            end
        end
    end

    return y
end
