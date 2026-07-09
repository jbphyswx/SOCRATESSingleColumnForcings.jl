#=
Conservative regridding: mass-matrix construction, spline-value inversion, non-negative least
squares, and the default conservative-interpolation kwargs.

Depends on `interpolants.jl` (included before this file) for `build_spline`, `safe_integrate`,
`FastLinear1DInterpolation`, and `AbstractInterpolationMethod`, and on `boundary_conditions.jl`
for the boundary-condition types. `LinearAlgebra`, `Statistics`, and `StaticArrays` are provided
by the enclosing `Interpolation` module.
=#

# Submodule-local NaN scrub (kept independent of the parent package so the submodule has no
# upward dependency). Used to clean leaky inversions in `conservative_spline_values`.
resolve_nan(x::FT, val::FT = zero(FT)) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0
function resolve_nan!(x::AbstractArray{FT}, val::FT = FT(0.0)) where {FT}
    @inbounds for i in eachindex(x)
        x[i] = resolve_nan(x[i], val)
    end
end

# non-negative least squares (used for enforce_positivity). Only a stub here; the sole method is
# provided by the NonNegLeastSquares extension. Calling this with `NonNegLeastSquares` unloaded raises
# a MethodError whose hint points at the extension's trigger package.
function nnls_solve end

#=
For the nnls_alg,
[https://www.sciencedirect.com/science/article/pii/S1877050917307858] explains the active and passive sets well.
[https://conservancy.umn.edu/server/api/core/bitstreams/ef5b0697-b7ea-4e5b-9d97-92b29a2f468e/content] chapter 2 points out that the pivot methods tend to be computationally optimal for NNLS (Kim and Park is in NonNegLeastSquares.jl)


The true profile shouldn't be incredibly spiky , and thus A\b should not induce many negative numbers in the A\b inversion). Thus the `active set` (numbers set to 0 because they violated the nonnegativity constraint) is usually small, so we usually don't need many iterations to find a solution

Because A is usually tridiagonal (maybe septadiagonal at worse for cubic splines) the inversions which might ordinarily be expensive are actually quite fast -- caching doesn't seem to improve performance.


Some sample timings with A = (1000x1000) and tridiagonal (but full matrix type, not LinearAlgebra.Tridiagonal) and M = (1000,) is fully dense and not sparse (few zeros):


@btime nonneg_lsq($A, $M; alg=:nnls) evals=2 samples=2  # doesn't support `tol`
    1.027 s (13 allocations: 7.68 MiB)

@btime nonneg_lsq($A, $M; alg=:nnls) evals=2 samples=2 


@btime nonneg_lsq($A, $M; alg=:fnnls, tol=1e-12) evals=2 samples=2 
    6.200 s (22914 allocations: 3.71 GiB)
@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-8) evals=2 samples=2
    6.258 s (22914 allocations: 3.71 GiB)


@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-8) evals=2 samples=2
    33.881 ms (48 allocations: 38.25 MiB)     : This is 
@btime nonneg_lsq($A, $M; alg=:nnls, tol=1e-12) evals=2 samples=2
    35.031 ms (48 allocations: 38.25 MiB)

So, :pivot is over 150x faster, with about 100x fewer allocations than :fnnls, but does allocate more than :mmls

For the pivot variants, :cache and :comb

@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-12, variant=:cache) evals=2 samples=2
    909.911 ms (67 allocations: 76.54 MiB)

@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-12, variant=:comb) evals=2 samples=2
    932.510 ms (113 allocations: 76.54 MiB)

I thk the active sets are just too small for caching things to be worth it


The only quirk is that :pivot can leave some very small negative numbers in the solution, which we have to resolve to 0, but that's very fast to do.

This probably also explains why :pivot was the default in the original code.

Note if you're doing fully random matrices for A, :pivot was slower, but then that's not really a realistic representation of what the workflow would be for regridding input variables like this where A is sparse and has a real structure from the grid.
=#

"""
    default_conservative_interp_kwargs

Default keyword bundle for conservative vertical regridding (`preserve_monotonicity`,
`enforce_positivity`, `nnls_alg`, etc.). Merge with overrides via
[`get_conservative_interp_kwargs`](@ref).
"""
const default_conservative_interp_kwargs = (;
    preserve_monotonicity = true,
    enforce_positivity = false,
    nnls_alg = :pivot, # should be optimal for most cases, see the discussion above.
    nnls_tol = 1e-8, # default for package
    enforce_conservation = true,
    integrate_method = :invert,
)

const DCIKT = typeof(default_conservative_interp_kwargs)
const DCIKDT = Dict{Symbol, Union{Bool, Symbol, Float64}} # Dict for conservative interpolation kwargs
const default_conservative_interp_kwargs_dict = DCIKDT(pairs(default_conservative_interp_kwargs))

get_conservative_interp_kwargs(::Nothing) = default_conservative_interp_kwargs

get_conservative_interp_kwargs(conservative_interp_kwargs::NamedTuple) = NamedTuple((
    key => get(conservative_interp_kwargs, key, default_conservative_interp_kwargs[key]) for
    key in keys(default_conservative_interp_kwargs)
)) # convert arbitrary NamedTuple or Dict to DCIKT, keeping only relevant keys and defaulting to the new tuple but falling back to the old one if not found

get_conservative_interp_kwargs(conservative_interp_kwargs::Dict) = get_conservative_interp_kwargs(
    NamedTuple((Symbol(k) => isa(v, String) ? Symbol(v) : v for (k, v) in conservative_interp_kwargs)),
) # convert any strings to symbols and then passed to the NamedTuple constructor


get_conservative_interp_kwargs(conservative_interp_kwargs::DCIKT) =
    merge(default_conservative_interp_kwargs, conservative_interp_kwargs) # merge between two DCIKTs (should maybe be faster than the iterative method above)

"""
    get_conservative_interp_kwargs(kwargs...)

Normalize user conservative-regrid kwargs against [`default_conservative_interp_kwargs`](@ref).
"""
get_conservative_interp_kwargs(; kwargs...) =
    get_conservative_interp_kwargs(merge(default_conservative_interp_kwargs, kwargs)) # merge kwargs (pairs iterator) into the default DCIKT then strip it down to the relevant keys / reorder if we added any


"""
    conservative_regridder(x, xp, yp; kwargs...)

Mass-conserving 1D regrid from source nodes `(xp, yp)` onto target coordinates `x`.
"""
function conservative_regridder(
    x::AbstractVector{FT},
    xp::AbstractVector,
    yp::AbstractVector;
    bc::BCT = ExtrapolateBoundaryCondition(), # must be extrapolate to integrate outside the knots (boundary of the data)... Dierckx integral ignored silently, we use dierckx_safe_integrate() instead
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
    f_enhancement_factor = 1,
    f_p_enhancement_factor = 1,
    integrate_method::Symbol = :invert,
    rtol::FT = FT(1e-6),
    preserve_monotonicity::Bool = false,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::FT = FT(1e-8),
    enforce_conservation::Bool = true,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
) where {BCT <: ValidBoundaryConditions, FT}

    # Build the spline
    # each node gets the mean of the spline over its area of influence, which we'll define as being the nearest neighbor
    # at the lower and upper end we'll extrapolate as far as the previous edge

    # get the spline
    spl = build_spline(method, xp, yp; bc = bc)

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

        y[i] = safe_integrate(spl, x_edges[i], x_edges[i + 1]; bc = bc) / (x_edges[i + 1] - x_edges[i])
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
            y_mean = Statistics.mean(y)
            if (0 < y_mean < 1) # large means are fine; ignore small ones and skip when the mean is 0
                # values below ~2*eps(FT)^0.5 become unstable, so bring small means up to 1
                y /= y_mean # scale the cell mean to 1 (data with a very large dynamic range is not a good fit for this regridding)
            end
        end

        # solve for the cell values that reproduce the integrated masses `y`, in place into `y`
        # (allocation-free via the cached factorization `Af`)
        _, y = conservative_spline_values!(
            y,
            x_edges,
            y;
            k = k,
            method = method,
            bc = bc,
            enforce_positivity = enforce_positivity,
            nnls_alg = nnls_alg,
            nnls_tol = nnls_tol,
            A = A,
            Af = Af,
        )

        if enforce_positivity && (0 < y_mean < 1)
            y *= y_mean
        end

        if preserve_monotonicity # ensure our new values are within the bounds of their adjacent points in the original grid
            for i in 1:(length(x))


                # Find index of largest point in xp smaller than x_edges[i]
                i_low = searchsortedlast(xp, x_edges[i]) - (k - 1) # index of the largest xp point below x_edges[i]. There is already a span of 1 between i_low and i_high, so k-1 extra width suffices (at k=1 that still leaves ≥2 contributing points, since edges are unique)
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
            total = safe_integrate(build_spline(method, x, y; bc = bc), x_edges[1], x_edges[end]; bc = bc)
            if !iszero(total)
                y *= safe_integrate(spl, x_edges[1], x_edges[end]; bc = bc) / total # rescale to the original total mass
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

# In-place workhorse: solve `A · yc = mc` for the cell-center values, writing the result into the
# caller-provided `yc` (which may alias `mc`). The linear-solve path uses `ldiv!` with the cached
# factorization `Af`, so it is allocation-free; the positivity path uses NNLS (no in-place solver
# exists, so its broadcast still allocates the NNLS result). Returns `(xc, yc)`.
function conservative_spline_values!(
    yc::AbstractVector,
    xf::AbstractVector{FT},
    mc::AbstractVector;
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::Real = 1e-8,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing, # precomputed factorization of A (preferred): enables the allocation-free `ldiv!`
) where {FT <: Real, BCT <: ValidBoundaryConditions}
    xc = xc_from_xf(xf)
    if enforce_positivity
        isnothing(A) && (A = conservative_mass_matrix(xc; bc = bc, k = k, method = method)) # NNLS can't use a factorization
        if 0 < Statistics.mean(mc) < (2 * eps(FT)^0.5)
            @warn "mean(mc) = $(Statistics.mean(mc)) < [2 * eps($FT)^0.5 = $(2 * eps(FT)^0.5)]; this is very small, NNLS may arbitrarily converge to 0. Consider scaling your data to be larger."
        end
        yc .= max.(nnls_solve(A, mc; alg = nnls_alg, tol = nnls_tol)[:], zero(FT)) # bound by 0 (e.g. :pivot can leave 1e-17 underflow negatives)
    elseif Af isa LinearAlgebra.Factorization
        yc === mc ? LinearAlgebra.ldiv!(Af, yc) : LinearAlgebra.ldiv!(yc, Af, mc) # allocation-free solve into yc
    else
        Amat = A !== nothing ? A : (Af !== nothing ? Af : conservative_mass_matrix(xc; bc = bc, k = k, method = method))
        F = LinearAlgebra.factorize(Amat)
        yc === mc ? LinearAlgebra.ldiv!(F, yc) : LinearAlgebra.ldiv!(yc, F, mc)
    end

    if any(!isfinite, mc)
        error("Received invalid (non-finite) input in mc = $mc")
    elseif any(!isfinite, yc)
        resolve_nan!(yc) # set NaNs (from underflow-tiny inputs) to zero
    end
    return xc, yc
end

# Allocating convenience: allocate the cell-center output, then delegate to the in-place solve.
function conservative_spline_values(
    xf::AbstractVector,
    mc::AbstractVector{FT2};
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::Real = 1e-8,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing,
) where {FT2 <: Real, BCT <: ValidBoundaryConditions}
    yc = zeros(FT2, length(xf) - 1) # should use similar() here
    return conservative_spline_values!(yc, xf, mc; bc = bc, k = k, method = method,
        enforce_positivity = enforce_positivity, nnls_alg = nnls_alg, nnls_tol = nnls_tol, A = A, Af = Af)
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

# 2D (time-batched) in-place workhorse: each column of `mc` is a cell-mean profile; solve into the
# matching column of the caller-provided `yc` (which may alias `mc`). Linear path is `ldiv!` with the
# cached factorization (allocation-free); positivity path solves NNLS per contiguous mass-bearing
# block of rows. Returns `(xc, yc)`.
function conservative_spline_values!(
    yc::AbstractMatrix,
    xf::AbstractVector{FT},
    mc::AbstractMatrix;
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::Real = 1e-8,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing,
) where {FT <: Real, BCT <: ValidBoundaryConditions}
    xc = xc_from_xf(xf)
    if enforce_positivity
        isnothing(A) && (A = conservative_mass_matrix(xc; bc = bc, k = k, method = method))
        if any(x -> (0 < x < 2 * eps(FT)^0.5), Statistics.mean(mc, dims = 1))
            @warn "mean(mc) very small (< 2·eps($FT)^0.5); NNLS may arbitrarily converge to 0. Consider scaling your data to be larger."
        end
        # only solve rows (cells) that carry mass — cheap for mostly-zero fields (e.g. condensate)
        for region in contiguous_true_ranges((!iszero).(mc); dim = 1)
            A_sub = @view A[region, region]
            mc_sub = @view mc[region, :]
            yc[region, :] .= max.(nnls_solve(A_sub, mc_sub; alg = nnls_alg, tol = nnls_tol), zero(FT)) # bound by 0 (:pivot underflow)
        end
    elseif Af isa LinearAlgebra.Factorization
        yc === mc ? LinearAlgebra.ldiv!(Af, yc) : LinearAlgebra.ldiv!(yc, Af, mc) # allocation-free batched solve
    else
        Amat = A !== nothing ? A : (Af !== nothing ? Af : conservative_mass_matrix(xc; bc = bc, k = k, method = method))
        F = LinearAlgebra.factorize(Amat)
        yc === mc ? LinearAlgebra.ldiv!(F, yc) : LinearAlgebra.ldiv!(yc, F, mc)
    end

    if any(!isfinite, mc)
        error("Received invalid (non-finite) input in mc. NaN inputs are not supported for conservative regridding (they break the matrix solve).")
    elseif any(!isfinite, yc)
        @error "NaN values in yc from conservative solve" xf bc k enforce_positivity nnls_alg
        resolve_nan!(yc) # set NaNs (from underflow-tiny inputs) to zero
    end
    return xc, yc
end

# Allocating convenience (2D): allocate the cell-center matrix, then delegate to the in-place solve.
function conservative_spline_values(
    xf::AbstractVector,
    mc::AbstractMatrix{FT};
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
    enforce_positivity::Bool = false,
    nnls_alg::Symbol = :pivot,
    nnls_tol::Real = 1e-8,
    A::Union{AbstractMatrix, Nothing} = nothing,
    Af::Union{AbstractMatrix, LinearAlgebra.Factorization, Nothing} = nothing,
) where {FT <: Real, BCT <: ValidBoundaryConditions}
    yc = zeros(FT, length(xf) - 1, size(mc, 2)) # should use similar() here
    return conservative_spline_values!(yc, xf, mc; bc = bc, k = k, method = method,
        enforce_positivity = enforce_positivity, nnls_alg = nnls_alg, nnls_tol = nnls_tol, A = A, Af = Af)
end

"""
    conservative_mass_matrix(xc; method = FastLinear1DInterpolation, bc = ExtrapolateBoundaryCondition(), k = 1)

Build the mass matrix used by conservative vertical regridding on cell centers `xc`.
"""
function conservative_mass_matrix(
    xc::AbstractVector{FT};
    bc::BCT = ExtrapolateBoundaryCondition(),
    k::Int = 1,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation,
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
        # allocates; a non-allocating similar()-based version would be better
        A = zeros(Float64, n, n) # need a float type: can't factorize/integrate/transform in integer arithmetic
    else
        A = zeros(FT, n, n)
    end

    @inbounds for j in 1:n
        ej = zeros(FT, n)
        ej[j] = FT(1.0)
        φj = build_spline(method, xc, ej; bc = bc)


        width = k + 1 # support width of each basis function (keeps this O(kn) not O(n^2))
        i_start = max(1, j - width)
        i_end = min(n, j + width)
        @inbounds for i in i_start:i_end
            A[i, j] = safe_integrate(φj, xf[i], xf[i + 1]; bc = bc) / (xf[i + 1] - xf[i])
        end
    end
    @. A = max(A, FT(0)) # clamp to non-negative: subtraction underflow can leave tiny (~1e-17) negatives
    return A
end
