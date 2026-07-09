#=
Interpolation method types, node-backing helpers, the fast piecewise-linear interpolant, spline
construction/evaluation, node pruning, and safe integration.

Boundary-condition types live in `boundary_conditions.jl` (included before this file).
Conservative regridding lives in `conservative.jl` (included after this file).
`StaticArrays`, `LinearAlgebra`, and `Statistics` are provided by the enclosing `Interpolation`
module.
=#

# --- interpolation method + interpolant supertypes ------------------------------------------ #
"""Supertype for 1D interpolation method singletons (e.g. [`FastLinear1DInterpolation`](@ref))."""
abstract type AbstractInterpolationMethod end

"""Supertype for built 1D interpolants (e.g. [`Fast1DLinearInterpolant`](@ref))."""
abstract type AbstractInterpolant end

# --- interval search ------------------------------------------------------------------------ #
# Interval index `i` bracketing the query between nodes `i`/`i+1`, clamped to `[1, N-1]`. There is NO
# layout phantom type and NO runtime uniformity test: the cost follows the coordinate's TYPE, because
# `searchsortedlast` is already O(1) on an `AbstractRange` (exact) and O(log N) on a `Vector`/`SVector`.
# So a caller who wants O(1) evaluation stores the (known-uniform) coordinate as a range. `asc` marks an
# ascending grid; descending grids are scanned in place without reversing storage.
@inline function _find_interval(xp, x, asc::Bool)
    N = length(xp)
    if asc
        i = searchsortedlast(xp, x)
    else
        i = 1
        @inbounds while i < N && xp[i] > x
            i += 1
        end
        i -= 1
    end
    return clamp(i, 1, N - 1)
end

# --- node/value storage coercion (two orthogonal knobs: backing × element type) -------------- #
create_svector(x::AbstractVector{FT}) where {FT} = StaticArrays.SVector{length(x), FT}(x) # SVector from an array (length looked up at compile time)

# `coerce_vector(backing, eltype, v)` coerces a working vector into a built interpolant's storage.
# `backing` and `eltype` are INDEPENDENT; either may be `nothing` for passthrough (keep the input's
# backing / element type). This is how the pipeline enforces requested storage per coordinate and value
# (fixing the ignored-eltype bug) while staying fully general:
#   coerce_vector(Vector, Float32, v)       -> Vector{Float32}
#   coerce_vector(SVector, nothing, v)      -> SVector{len, eltype(v)}
#   coerce_vector(StepRangeLen, nothing, v) -> range over v (uniform axis; O(1) lookup)
#   coerce_vector(nothing, Float32, v)      -> keep v's backing, coerce eltype -> Float32
#   coerce_vector(nothing, nothing, v)      -> v unchanged (full passthrough)
# Add a `coerce_vector(::Type{MyBacking}, eltype, v)` method to support a custom backing.
# A storage spec is a `Tuple{Backing, Eltype}` TYPE (two knobs, bundled); `Nothing` (the type) means
# passthrough for either knob. It's a TYPE (not a value tuple) so it participates in specialization —
# `coerce_vector` infers its result, keeping the built interpolant (and the pipeline's returned NamedTuple)
# concretely typed. Taking the spec first makes it `Base.Fix1`-friendly (`Base.Fix1(coerce_vector, spec)`),
# matching the original one-arg-fixed idiom.
@inline _resolve_eltype(::Type{Nothing}, v::AbstractVector) = Base.eltype(v)
@inline _resolve_eltype(::Type{FT}, ::AbstractVector) where {FT} = FT

"""
    coerce_vector(spec, v)

Coerce vector `v` into the backing and element type requested by `spec::Type{Tuple{Backing, Eltype}}`.
Either knob may be `Nothing` for passthrough on that axis. Used by [`build_spline`](@ref) and
[`get_column_forcing`](@ref) storage specs.
"""
coerce_vector(::Type{Tuple{B, E}}, v::AbstractVector) where {B, E} = _coerce_storage(B, E, v)

_coerce_storage(::Type{Nothing}, ::Type{E}, v::AbstractVector) where {E} =
    (T = _resolve_eltype(E, v); T === Base.eltype(v) ? v : T.(v))  # keep backing family (broadcast), coerce eltype
_coerce_storage(::Type{Vector}, ::Type{E}, v::AbstractVector) where {E} = convert(Vector{_resolve_eltype(E, v)}, v)
_coerce_storage(::Type{<:StaticArrays.SVector}, ::Type{E}, v::AbstractVector) where {E} =
    StaticArrays.SVector{length(v), _resolve_eltype(E, v)}(v)
_coerce_storage(::Type{<:AbstractRange}, ::Type{E}, v::AbstractVector) where {E} = _coerce_uniform_range(_resolve_eltype(E, v), v)

# Build an exactly-uniform range spanning `v` (for a known-uniform axis, e.g. the Atlas hourly time axis).
# Uses EXACT step equality (integer/exact), NOT a float tolerance — the range makes interval search O(1).

function _coerce_uniform_range(::Type{T}, v::AbstractVector{FT}, check_values::Bool = true; tol::FT = zero(FT)) where {T, FT} # could use ::Val{check_values} instead of check_values::Bool
    n = length(v)
    n <= 1 && error("coerce_vector to a range needs ≥ 2 nodes (got $n)")
    @inbounds x0 = v[1]
    @inbounds s = v[2] - v[1]
    if check_values
        @inbounds for i in 3:n
            if iszero(tol)
                (v[i] - v[i - 1]) == s || error("coerce_vector to a range requires an exactly-uniform axis (non-uniform spacing at index $i)")
            else
                abs(v[i] - v[i - 1]) - s <= tol || error("coerce_vector to a range requires an exactly-uniform axis (non-uniform spacing at index $i)")
            end
        end
    end
    return range(T(x0); step = T(s), length = n)
end
function _coerce_uniform_range(v::AbstractVector{FT}, check_values::Bool = true; tol::FT = zero(FT)) where {FT} # could use ::Val{check_values} instead of check_values::Bool
    return _coerce_uniform_range(eltype(v), v, check_values; tol = tol)
end

@inline _coerce_uniform_range(v::AbstractRange) = v
function _coerce_uniform_range(::Type{T}, v::AbstractRange{FT}, check_values::Bool = true; tol::FT = zero(FT)) where {T, FT} # could use ::Val{check_values} instead of check_values::Bool
    return range(T(v.start); step = T(v.step), length = length(v))
end

# --- fast linear interpolation method ------------------------------------------------------- #
"""Piecewise-linear 1D interpolation method (fast path)."""
struct FastLinear1DInterpolationMethod <: AbstractInterpolationMethod end

"""Default piecewise-linear interpolation method instance."""
const FastLinear1DInterpolation = FastLinear1DInterpolationMethod()

# Core scalar evaluation. `asc` is derived from the stored nodes so ascending and descending grids share
# one code path; the interval search's cost follows `typeof(xp)` (range → O(1)). A single node is constant
# (degenerate the eval formula so the return type stays stable across the N==1 / N>1 branches).
@inline function _eval_linear(xp::AbstractVector, fp::AbstractVector, bc, x)
    N = length(xp)
    @inbounds x0n = xp[1]
    if isone(N)
        y1 = convert(promote_type(eltype(fp), eltype(xp), typeof(x)), @inbounds fp[1]) # type sability
        if (bc isa NearestBoundaryCondition) || (bc isa ExtrapolateBoundaryCondition)
            return y1
        else
            return (x == x0n) ? (y1) : throw(BoundsError(xp, x))
        end
    end
    @inbounds xmin, xmax = xp[1], xp[N]
    asc = xmin <= xmax

    # --- Out-of-bounds ---
    if asc
        if x < xmin
            if bc isa NearestBoundaryCondition
                return @inbounds fp[1] + (fp[1] - fp[1]) * (x - x0n) / oneunit(x - x0n)
            elseif bc isa ExtrapolateBoundaryCondition
                return @inbounds fp[1] + (fp[2] - fp[1]) * (x - xp[1]) / (xp[2] - xp[1])
            else
                error("x = $x below interpolation range [$xmin, $xmax]")
            end
        elseif x > xmax
            if bc isa NearestBoundaryCondition
                return @inbounds fp[N] + (fp[N] - fp[N]) * (x - x0n) / oneunit(x - x0n)
            elseif bc isa ExtrapolateBoundaryCondition
                return @inbounds fp[N - 1] + (fp[N] - fp[N - 1]) * (x - xp[N - 1]) / (xp[N] - xp[N - 1])
            else
                error("x = $x above interpolation range [$xmin, $xmax]")
            end
        end
    else
        if x > xmin
            if bc isa NearestBoundaryCondition
                return @inbounds fp[1] + (fp[1] - fp[1]) * (x - x0n) / oneunit(x - x0n)
            elseif bc isa ExtrapolateBoundaryCondition
                return @inbounds fp[1] + (fp[2] - fp[1]) * (x - xp[1]) / (xp[2] - xp[1])
            else
                error("x = $x above interpolation range [$xmin, $xmax]")
            end
        elseif x < xmax
            if bc isa NearestBoundaryCondition
                return @inbounds fp[N] + (fp[N] - fp[N]) * (x - x0n) / oneunit(x - x0n)
            elseif bc isa ExtrapolateBoundaryCondition
                return @inbounds fp[N - 1] + (fp[N] - fp[N - 1]) * (x - xp[N - 1]) / (xp[N] - xp[N - 1])
            else
                error("x = $x below interpolation range [$xmin, $xmax]")
            end
        end
    end

    # --- Interval search (cost follows typeof(xp)) + linear interpolation ---
    i = _find_interval(xp, x, asc)
    @inbounds x0, x1 = xp[i], xp[i + 1]
    @inbounds y0, y1 = fp[i], fp[i + 1]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

# Scalar fast linear interpolation (free function on raw node/value vectors).
function fast1d_interp(
    x::TT,
    xp::AbstractVector{T},
    fp::AbstractVector;
    bc::BCT = ErrorBoundaryCondition(),
    monotonic::Bool = false,
    check_monotonic::Bool = false,
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}
    if monotonic && check_monotonic
        N = length(xp)
        asc = true
        desc = true
        @inbounds for i in 2:N
            xp[i] < xp[i - 1] && (asc = false)
            xp[i] > xp[i - 1] && (desc = false)
            (!asc && !desc) && error("xp must be monotonic (ascending or descending)")
        end
    end
    return _eval_linear(xp, fp, bc, x)
end

function fast1d_interp!(
    y::AbstractVector{TT},
    x::AbstractVector{TT},
    xp::AbstractVector{T},
    fp::AbstractVector;
    bc::BCT = ErrorBoundaryCondition(),
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}
    @inbounds for j in eachindex(x)
        y[j] = _eval_linear(xp, fp, bc, x[j])
    end
    return y
end

# Vector version (broadcastable)
function fast1d_interp(
    x::AbstractVector{TT},
    xp::AbstractVector{T},
    fp::AbstractVector;
    bc::BCT = ErrorBoundaryCondition(),
) where {T <: Real, TT <: Real, BCT <: ValidBoundaryConditions}
    y = similar(x, promote_type(TT, T, eltype(fp)))
    fast1d_interp!(y, x, xp, fp; bc = bc)
    return y
end

# --- Fast1DLinearInterpolant ---------------------------------------------------------------- #
# Nodes (`xp::X`) and values (`fp::F`) have INDEPENDENT backings and element types (decoupled). The node
# layout is encoded by `X` itself: a range coordinate ⇒ O(1) evaluation, a `Vector`/`SVector` ⇒ O(log N).
struct Fast1DLinearInterpolant{X <: AbstractVector, F <: AbstractVector, BCT <: ValidBoundaryConditions} <: AbstractInterpolant
    xp::X
    fp::F
    bc::BCT
end

# `drop_collinear` is a `Val` so the choice is a compile-time constant: `Val(false)` dispatches to the
# identity (no pruning) — type-stable, allocation-free, and required for a range `xp`. `Val(true)` is the
# cold pruning path (its result length/type depends on the data).
@inline _maybe_prune(::Val{false}, xp, fp) = (xp, fp)
@inline _maybe_prune(::Val{true}, xp, fp) = length(xp) > 2 ? drop_collinear_nodes(xp, fp) : (xp, fp)

"""
    Fast1DLinearInterpolant(xp, fp; bc, drop_collinear = Val(false))

Build a piecewise-linear interpolant on nodes `(xp, fp)` (equal length). This primitive is
backing/eltype-PRESERVING — it stores exactly the `xp`/`fp` it is given (a range `xp` stays a range,
giving O(1) evaluation). Choosing storage (backing + element type per coord/value) is the caller's job,
done via [`coerce_vector`](@ref) at the pipeline layer.

When `drop_collinear = Val(true)`, collinear interior nodes are pruned via [`drop_collinear_nodes`](@ref), which
preserves each node backing (`Vector -> Vector`, `SVector{N} -> SVector{keep}`). The pruned length depends
on the node *data*, so for length-typed backings (`SVector`) construction is type-unstable and allocates
at build time (a cold path). Use `drop_collinear = Val(false)` for the non-allocating, type-stable path (and
required for a range `xp`, since pruning would break the range).
"""
function Fast1DLinearInterpolant(
    xp::AbstractVector,
    fp::AbstractVector;
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Val = Val(false),
) where {BCT <: ValidBoundaryConditions}
    @assert length(xp) == length(fp) "xp and fp must have the same length"
    xk, fk = _maybe_prune(drop_collinear, xp, fp)
    return Fast1DLinearInterpolant{typeof(xk), typeof(fk), BCT}(xk, fk, bc)
end

# Make it broadcastable (for vector or static array evaluation)
Base.broadcastable(s::Fast1DLinearInterpolant) = tuple(s)

# ---- pretty printing ---------------------------------------------------------------------------- #
function Base.show(io::IO, itp::Fast1DLinearInterpolant)
    n = length(itp.xp)
    print(io, "Fast1DLinearInterpolant(", n, " node", n == 1 ? "" : "s")
    n > 0 && print(io, ", x∈[", first(itp.xp), ", ", last(itp.xp), "]")
    print(io, ", ", bc_string(itp.bc), " bc)")
end

function Base.show(io::IO, ::MIME"text/plain", itp::Fast1DLinearInterpolant)
    n = length(itp.xp)
    println(io, "Fast1DLinearInterpolant with ", n, " node", n == 1 ? "" : "s", ":")
    ioc = IOContext(io, :limit => true, :compact => true)
    print(ioc, "  xp: ", length(itp.xp), "-element ", typeof(itp.xp), ": "); show(ioc, itp.xp); println(ioc)
    print(ioc, "  fp: ", length(itp.fp), "-element ", typeof(itp.fp), ": "); show(ioc, itp.fp); println(ioc)
    print(ioc, "  bc: ", bc_string(itp.bc))
end
# --------------------------------------------------------------------------------------------------- #

change_bc(
    sp::Fast1DLinearInterpolant{X, F, BCTO},
    new_bc::BCTN,
) where {X, F, BCTO <: ValidBoundaryConditions, BCTN <: ValidBoundaryConditions} =
    Fast1DLinearInterpolant{X, F, BCTN}(sp.xp, sp.fp, new_bc)

# Callable
(s::Fast1DLinearInterpolant{T,S,B})(x) where {T,S,B} = _eval_linear(s.xp, s.fp, s.bc, x)
function fast1d_derivative(
    x::T,
    spl::Fast1DLinearInterpolant;
) where {T <: Real}
    xp, fp = spl.xp, spl.fp
    N = length(xp)
    bc = spl.bc
    @inbounds asc = xp[1] <= xp[N]

    below = asc ? x < xp[1] : x > xp[1]
    above = asc ? x > xp[N] : x < xp[N]
    if below
        if bc isa NearestBoundaryCondition
            return zero(eltype(fp)) / oneunit(eltype(xp))
        elseif bc isa ExtrapolateBoundaryCondition
            return @inbounds (fp[2] - fp[1]) / (xp[2] - xp[1])
        else
            error("x = $x outside interpolation range and bc=$bc")
        end
    elseif above
        if bc isa NearestBoundaryCondition
            return zero(eltype(fp)) / oneunit(eltype(xp))
        elseif bc isa ExtrapolateBoundaryCondition
            return @inbounds (fp[N] - fp[N - 1]) / (xp[N] - xp[N - 1])
        else
            error("x = $x outside interpolation range and bc=$bc")
        end
    else
        i = _find_interval(xp, x, asc)
        return @inbounds (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
    end
end

function fast1d_derivative!(
    y::AbstractVector, # output array
    x::AbstractVector, # derivative locations
    spl::Fast1DLinearInterpolant,
)
    @inbounds for j in eachindex(x)
        y[j] = fast1d_derivative(x[j], spl)
    end
    return y
end
# Vectorized derivative
function fast1d_derivative(
    x::AbstractVector{T},
    spl::Fast1DLinearInterpolant,
) where {T <: Real}
    y = similar(x, promote_type(T, eltype(spl.fp)))
    fast1d_derivative!(y, x, spl)
    return y
end


# ============================================================================================================================================= #
# Custom Types
# ============================================================================================================================================= #

"""
    UniformRange(start, step, length)

A custom range type that is exactly uniform.
Julia default StepRangeLen comes with twiceprecision arithmtic that massively inflates inference times.
TODO :: Agent fill in the rest of this and the actual implementation
"""
struct UniformRange{T, LT <: Integer} <: AbstractRange{T}
    start::T
    step::T
    length::LT
    inv_step::T
end
UniformRange(start::T, step::T, length::LT) where {T, LT <: Integer} = UniformRange{T, LT}(start, step, length, inv(step))
Base.length(r::UniformRange) = Int(r.length)
Base.size(r::UniformRange) = (Int(r.length),)
Base.first(r::UniformRange) = r.start
Base.step(r::UniformRange)  = r.step
Base.last(r::UniformRange)  = r.start + (r.length - 1) * r.step
@inline function Base.getindex(r::UniformRange, i::Int)
    @boundscheck (1 <= i <= r.length) || throw(BoundsError(r, i))
    return r.start + (i - 1) * r.step
end
@inline function Base.getindex(r::UniformRange, s::AbstractRange{<:Integer})
    @boundscheck checkbounds(r, s)
    return UniformRange(r.start + (first(s) - 1) * r.step, r.step * step(s), length(s))
end


function _eval_linear(xp::UniformRange, fp::AbstractVector, bc::BCT, x) where {BCT <: ValidBoundaryConditions}
    N = xp.length
    x0 = xp.start
    h  = xp.step
    if isone(N)
        if (bc isa NearestBoundaryCondition) || (bc isa ExtrapolateBoundaryCondition)
            return @inbounds fp[1]
        else
            return (x == x0) ? (@inbounds fp[1]) : throw(BoundsError(xp, x))
        end
    end
    xmin, xmax = x0, x0 + (N - 1) * h
    asc = xmin <= xmax

    if asc
        if x < xmin
            bc isa NearestBoundaryCondition     && return @inbounds fp[1]
            bc isa ExtrapolateBoundaryCondition && return @inbounds fp[1] + (fp[2] - fp[1]) * (x - xmin) * xp.inv_step
            error("x = $x below interpolation range [$xmin, $xmax]")
        elseif x > xmax
            bc isa NearestBoundaryCondition     && return @inbounds fp[N]
            bc isa ExtrapolateBoundaryCondition && return @inbounds fp[N - 1] + (fp[N] - fp[N - 1]) * (x - (xmax - h)) * xp.inv_step
            error("x = $x above interpolation range [$xmin, $xmax]")
        end
    else
        if x > xmin
            bc isa NearestBoundaryCondition     && return @inbounds fp[1]
            bc isa ExtrapolateBoundaryCondition && return @inbounds fp[1] + (fp[2] - fp[1]) * (x - xmin) * xp.inv_step
            error("x = $x above interpolation range [$xmax, $xmin]")
        elseif x < xmax
            bc isa NearestBoundaryCondition     && return @inbounds fp[N]
            bc isa ExtrapolateBoundaryCondition && return @inbounds fp[N - 1] + (fp[N] - fp[N - 1]) * (x - (xmax - h)) * xp.inv_step
            error("x = $x below interpolation range [$xmax, $xmin]")
        end
    end

    # in-bounds: pure multiply, index only fp (works for ascending and descending since t >= 0 in range)
    t = (x - x0) * xp.inv_step
    i = clamp(unsafe_trunc(Int, t) + 1, 1, N - 1)
    @inbounds y0 = fp[i]
    @inbounds y1 = fp[i + 1]
    return y0 + (y1 - y0) * (t - (i - 1))
end

function _coerce_storage(::Type{UniformRange}, ::Type{E}, v::AbstractVector) where {E}
    T = _resolve_eltype(E, v)
    n = length(v)
    n >= 2 || error("UniformRange needs ≥ 2 nodes (got $n)")
    @inbounds x0 = v[1]
    @inbounds s = v[2] - v[1]
    @inbounds for i in 3:n
        (v[i] - v[i-1]) == s || error("UniformRange requires an exactly-uniform axis (non-uniform at index $i)")
    end
    return UniformRange(T(x0), T(s), n)

end

# =============================== #
# =============================== #

"""Single-value `AbstractVector` backing for exactly-constant fields."""
struct Constant{T} <: AbstractVector{T}  # length-1 backing (elide boundschecking)
    value::T
end
Base.IndexStyle(::Type{<:Constant}) = IndexLinear()
Base.size(::Constant{T}) where {T} = (1,)
Base.axes(::Constant{T}) where {T} = (Base.OneTo(1),)
Base.length(::Constant{T}) where {T} = 1
@inline Base.getindex(c::Constant, ::Int) = c.value
(itp::Fast1DLinearInterpolant{T, <: Constant, BCTO} where {T, BCTO <: ValidBoundaryConditions})(x) = begin
    if (itp.bc isa NearestBoundaryCondition) || (itp.bc isa ExtrapolateBoundaryCondition)
        return itp.fp.value
    else
        return (x == itp.xp.value) ? itp.fp.value : throw(BoundsError(itp.xp, x))
    end
end

function _coerce_storage(::Type{Constant}, ::Type{E}, v::AbstractVector) where {E}
    T = _resolve_eltype(E, v)
    @inbounds all(==(v[1]), v) || error("Constant backing requires all-equal values")
    return Constant{T}(T(v[1]))
end

"""
    ConstantVector(value)

Length-`N` backing where every element equals `value` (for exactly-constant fields).
"""
struct ConstantVector{T, N} <: AbstractVector{T}
    value::T
end
ConstantVector(value::T) where {T} = ConstantVector{T, 1}(value)
ConstantVector(value::T, N::Integer) where {T} = ConstantVector{T, N}(value)
Base.IndexStyle(::Type{<:ConstantVector}) = IndexLinear()
Base.size(::ConstantVector{T, N}) where {T, N} = (N,)
Base.axes(::ConstantVector{T, N}) where {T, N} = (Base.OneTo(N),)
Base.length(::ConstantVector{T, N}) where {T, N} = N
@inline Base.getindex(c::ConstantVector, i::Int) = (firstindex(c) ≤ i ≤ lastindex(c)) ? c.value : throw(BoundsError(c, i))
function _coerce_storage(::Type{ConstantVector}, ::Type{E}, v::AbstractVector) where {E}
    T = _resolve_eltype(E, v)
    @inbounds all(==(v[1]), v) || error("Constant backing requires all-equal values")
    return ConstantVector{T, length(v)}(T(v[1]))
end
(itp::Fast1DLinearInterpolant{T, <: ConstantVector, BCTO} where {T, BCTO <: ValidBoundaryConditions})(x) = begin
    if (itp.bc isa NearestBoundaryCondition) || (itp.bc isa ExtrapolateBoundaryCondition)
        return itp.fp.value
    else
        return (itp.xp[begin] <= x <= itp.xp[end]) ? itp.fp.value : throw(BoundsError(itp.xp, x))
    end
end
# ============================================================================================================================================= #

function ConstantVector(value::Fast1DLinearInterpolant, ::Val{perform_checks} = Val(true)) where {perform_checks}
    if perform_checks
        all(isequal(value.fp[1]), value.fp) || error("ConstantVector requires all-equal values")
    end
    return ConstantVector(value.fp[1], length(value.xp))
end

function Constant(value::Fast1DLinearInterpolant, ::Val{perform_checks} = Val(true)) where {perform_checks}
    if perform_checks
        all(isequal(value.fp[1]), value.fp) || error("Constant requires all-equal values")
    end
    return Constant(value.fp[1])
end


function constantize_interpolant(itp::Fast1DLinearInterpolant, new_bc::ValidBoundaryConditions = itp.bc)
    
    if (new_bc isa ExtrapolateBoundaryCondition) || (new_bc isa NearestBoundaryCondition)
        new_xp = Constant(itp.xp[1])
        new_fp = Constant(itp.fp[1])
    elseif new_bc isa ErrorBoundaryCondition
        new_xp = itp.xp
        new_fp = ConstantVector(itp.fp[1], length(itp.xp))
    else
        error("Unsupported boundary condition: $new_bc")
    end

    return Fast1DLinearInterpolant(new_xp, new_fp; bc = new_bc, drop_collinear = Val(false))
end

# ============================================================================================================================================= #

# Build the interpolant for `FastLinear1DInterpolationMethod`. Backing/eltype-preserving (stores what it's
# given); the caller controls storage via [`coerce_vector`](@ref) upstream.
"""
    build_spline(method, xp, fp; bc = ErrorBoundaryCondition(), drop_collinear = Val(true))

Build a 1D interpolant for `method` on nodes `(xp, fp)`. Default method is
[`FastLinear1DInterpolation`](@ref). `drop_collinear` is a `Val`, so `Val(false)` keeps construction
type-stable and allocation-free (and is required for a range-backed `xp`).
"""
function build_spline(
    ::FastLinear1DInterpolationMethod,
    xp::AbstractVector,
    fp::AbstractVector;
    bc::BCT = ErrorBoundaryCondition(),
    drop_collinear::Val = Val(true),
) where {BCT <: ValidBoundaryConditions}
    return Fast1DLinearInterpolant(xp, fp; bc = bc, drop_collinear = drop_collinear)
end

# --------------------------------------------------------------------------------------------------------------------------------------------- #

# Build the interpolant for `method` and evaluate it at `x`.
"""
    interpolate_1d(x, xp, fp, method = FastLinear1DInterpolation; bc = ErrorBoundaryCondition())

Build a spline from `(xp, fp)` and evaluate at `x` (scalar or array).
"""
function interpolate_1d(
    x,
    xp,
    fp,
    method::AbstractInterpolationMethod = FastLinear1DInterpolation;
    bc::BCT = ErrorBoundaryCondition(),
) where {BCT <: ValidBoundaryConditions}
    if (xp isa AbstractVector) && (fp isa AbstractVector)
        xpT = Base.nonmissingtype(eltype(xp))
        fpT = Base.nonmissingtype(eltype(fp))
        if (xpT != eltype(xp)) || (fpT != eltype(fp))
            valid = .!ismissing.(xp) .& .!ismissing.(fp)
            xp = xp[valid]
            fp = fp[valid]
            T = promote_type(xpT, fpT)
            xp = T.(xp)
            fp = T.(fp)
        end
    end
    spl = build_spline(method, xp, fp; bc = bc)
    return spl.(x)
end


# ============================================================================================================================================= #
# Node pruning / shared-node coercion (operate on the interpolant's stored nodes)

# exact, division-free collinearity test for the middle point (b) of (a,b,c):
# slope(a,b) == slope(b,c)  <=>  (yb-ya)*(xc-xb) == (yc-yb)*(xb-xa)
@inline _is_collinear(xa, ya, xb, yb, xc, yc) = (yb - ya) * (xc - xb) == (yc - yb) * (xb - xa)

# Keep node `i` of `n`? Endpoints always; interior nodes only if they bend the line. Shared by every
# `drop_collinear_nodes` method (the `i==1 || i==n` short-circuit guards the neighbor accesses).
@inline _keep_node(x, y, i, n) =
    i == 1 || i == n || !_is_collinear(x[i - 1], y[i - 1], x[i], y[i], x[i + 1], y[i + 1])

"""
    drop_collinear_nodes(x, y) -> (x_kept, y_kept)

Drop interior nodes that lie exactly on the line through their neighbors (redundant for a
piecewise-linear interpolant). Returns the inputs unchanged when nothing is collinear (`n <= 2` or a
fully-bent polyline).

Backing-preserving PER ARRAY: `x` and `y` are pruned independently, each staying in its own backing
family — `Vector -> Vector`, `SVector{N} -> SVector{keep}` (isbits kept) — so one side's backing never
forces the other's (fixing the old coupling that demoted both to `Vector` for any mismatched pair). The
collinearity keep-mask depends on BOTH `x` and `y`, so it is threaded into each per-array prune. A range
(`AbstractRange`) coordinate cannot represent an irregular pruned node set: it is rebuilt via
`_coerce_uniform_range`, which returns a range when the kept nodes remain exactly uniform (e.g.
only the endpoints survive) and errors otherwise — the honest boundary instead of a silent demotion to
an allocating `Vector`.
"""
function drop_collinear_nodes(x::AbstractVector, y::AbstractVector)
    n = length(x)
    @assert length(y) == n "x and y must have equal length"
    n <= 2 && return (x, y)

    # count kept nodes once (endpoints + interior bends); this mask is shared by both per-array prunes
    keep = 0
    @inbounds for i in 1:n
        _keep_node(x, y, i, n) && (keep += 1)
    end
    keep == n && return (x, y) # nothing collinear -> both backings preserved trivially

    # materialize each side INDEPENDENTLY in its own backing (x's backing never forces y's, and vice versa)
    return (_prune_backing(x, x, y, keep), _prune_backing(y, x, y, keep))
end

# Prune array `v` to the `keep` nodes kept by the `(x, y)` collinearity mask, dispatched on `v`'s backing
# so `x` and `y` are pruned independently in their own storage families.

# default (Vector, and any backing whose `similar(v, keep)` returns its own type): exact-size fill.
function _prune_backing(v::AbstractVector, x, y, keep)
    n = length(v)
    out = similar(v, keep)
    @inbounds begin
        j = 0
        for i in 1:n
            if _keep_node(x, y, i, n)
                j += 1
                out[j] = v[i]
            end
        end
    end
    return out
end

# SVector: prune to `SVector{keep}` (keeping the isbits static backing). `keep` is data-dependent, so
# `SVector{keep}` is a runtime-typed result — this build is type-unstable and allocates at build time (a
# cold path), but the payoff is a SHORTER static vector: cheaper at use time and kinder to the compiler
# than a long one. The `Val(keep)` dispatch is the single, unavoidable runtime->type boundary; `sacollect`
# then gathers the kept indices into `SVector{K,Int}` and SVector-indexing returns `SVector{K}`.
_prune_backing(v::StaticArrays.SVector, x, y, keep) = _prune_svector(Val(keep), v, x, y)
@inline function _prune_svector(::Val{K}, v, x, y) where {K}
    n = length(v)
    idx = StaticArrays.sacollect(StaticArrays.SVector{K, Int}, i for i in 1:n if _keep_node(x, y, i, n))
    return v[idx]
end

# Range: an irregular pruned node set is not representable as a range. `v` is uniform, so uniform VALUE
# spacing ⟺ uniform INDEX-gap spacing of the kept nodes; check the kept-index gaps directly and return the
# sub-range `v[i0:Δ:prev]` (still a range, NO allocation), erroring on a non-uniform gap instead of silently
# demoting to a `Vector`. Endpoints are always kept, so `keep >= 2` in practice (`keep` is unused here — the
# gaps are re-derived in one pass — but the signature matches the other `_prune_backing` methods). The
# empty/singleton returns are defensive.
function _prune_backing(v::AbstractRange, x, y, keep)
    n = length(v)
    i0 = 0    # first kept index
    prev = 0  # previous kept index
    Δ = 0     # constant kept-index gap
    @inbounds for i in 1:n
        _keep_node(x, y, i, n) || continue
        if i0 == 0
            i0 = prev = i
        elseif Δ == 0
            Δ = i - prev
            prev = i
        else
            i - prev == Δ ||
                error("drop_collinear_nodes: pruned range nodes are not uniformly spaced (cannot stay an AbstractRange); use a Vector/SVector backing or drop_collinear = Val(false)")
            prev = i
        end
    end
    i0 == 0 && return v[1:0]     # empty (defensive; endpoints are always kept)
    Δ == 0 && return v[i0:i0]    # singleton (defensive)
    return v[i0:Δ:prev]
end

function drop_collinear_nodes(::AbstractInterpolant)
    error("Not implemented")
end

# Prune collinear interior nodes of a built linear interpolant, returning an equivalent interpolant
# on the reduced node set (nodes already pruned, so no second pass).
function drop_collinear_nodes(itp::Fast1DLinearInterpolant)
    xp_kept, fp_kept = drop_collinear_nodes(itp.xp, itp.fp)
    return Fast1DLinearInterpolant(xp_kept, fp_kept; bc = itp.bc, drop_collinear = Val(false))
end


"""
    coerce_to_shared_nodes(itp_collection)

Re-evaluate a collection of [`Fast1DLinearInterpolant`](@ref)s on a shared sorted node set.
"""
function coerce_to_shared_nodes(itp_collection)
    error("Not implemented")
end

# Rebuild every linear interpolant in the collection on the shared, sorted union of all their
# nodes. Each rebuilt interpolant keeps its own boundary condition and retains every shared node
# (`drop_collinear = Val(false)`) so the collection is defined on one common node vector.
function coerce_to_shared_nodes(itp_collection::AbstractVector{T}) where {T <: Fast1DLinearInterpolant}
    xs = sort(unique(reduce(vcat, (itp.xp for itp in itp_collection))))
    return [Fast1DLinearInterpolant(xs, itp.(xs); bc = itp.bc, drop_collinear = Val(false)) for itp in itp_collection]
end

function coerce_to_shared_nodes(itp_collection::StaticArrays.SVector{N, T}) where {N, T <: Fast1DLinearInterpolant}
    xs = sort(unique(reduce(vcat, (itp.xp for itp in itp_collection))))
    return StaticArrays.SVector{N}(
        Fast1DLinearInterpolant(xs, itp.(xs); bc = itp.bc, drop_collinear = Val(false)) for itp in itp_collection
    )
end

function coerce_to_shared_nodes(itp_collection::NTuple{N, <:Fast1DLinearInterpolant}) where {N}
    xs = sort(unique(reduce(vcat, (itp.xp for itp in itp_collection))))
    return ntuple(
        i -> Fast1DLinearInterpolant(xs, itp_collection[i].(xs); bc = itp_collection[i].bc, drop_collinear = Val(false)),
        Val(N),
    )
end


# ============================================================================================================================================= #
# Safe integration of the fast linear interpolant (bc respected directly in integration)

"""
Safe integration of Fast1DLinearInterpolant, respecting bc modes.

Unlike Dierckx, Fast1D respects boundary conditions in integration, so we do not pass bc in
"""
function fast1d_safe_integrate(
    spl::Fast1DLinearInterpolant,
    x1::T1,
    x2::T2;
) where {T1 <: Real, T2 <: Real}
    ST = eltype(spl.xp)
    FT = promote_type(T1, T2, ST, eltype(spl.fp))

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
    spl::Fast1DLinearInterpolant,
    x1::T1,
    x2::T2,
) where {T1 <: Real, T2 <: Real}
    ST = eltype(spl.xp)
    FT = promote_type(T1, T2, ST, eltype(spl.fp))

    xp, fp = spl.xp, spl.fp
    N = length(xp)
    y = zero(FT)

    # Loop over each interval intersecting [x1,x2] (grid assumed ascending here; cost follows typeof(xp))
    i1 = clamp(_find_interval(xp, x1, true), 1, N - 1)
    i2 = clamp(_find_interval(xp, x2, true), 1, N - 1)

    for i in i1:i2
        xa = max(x1, xp[i])
        xb = min(x2, xp[i + 1])
        slope = (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
        # ∫ (fp[i] + slope·(x - xp[i])) dx over [xa, xb]. The quadratic term is measured from the
        # segment's left node xp[i], NOT from xa — using (xb-xa)^2 would only be correct when xa == xp[i]
        # and under-counts the lower partial cell when x1 falls inside a segment.
        y += fp[i] * (xb - xa) + slope / 2 * ((xb - xp[i])^2 - (xa - xp[i])^2)
    end

    return y
end

"""
BC is only passed for Dierckx splines, Fast1D respects bc in integration already.
"""
safe_integrate(spl, x1::FT, x2::FT; bc = ExtrapolateBoundaryCondition()) where {FT} =
    error("safe_integrate not implemented for spline type $(typeof(spl))")
safe_integrate(
    spl::Fast1DLinearInterpolant,
    x1::FT,
    x2::FT;
    bc = ExtrapolateBoundaryCondition(),
) where {FT} = fast1d_safe_integrate(spl, x1, x2) # We do not pass bc through
