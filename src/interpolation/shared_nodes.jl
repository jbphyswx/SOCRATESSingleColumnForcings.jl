# ============================================================================
# Sorted union of interpolation nodes
#
# Each input node collection must already be sorted.
#
# Binary-heap k-way merge:
#   time:   O(M log N)
#   memory: O(N + U)
#
# where
#   N = number of interpolants
#   M = total number of input nodes
#   U = number of unique output nodes
#
# Critically, this NEVER constructs the O(M) concatenation.
# ============================================================================

# Fast1DLinearInterpolant already exposes its nodes directly.
@inline _shared_node_source(itp::Fast1DLinearInterpolant) = itp.xp

# Generic interpolants use the public/general node accessor.
@inline _shared_node_source(itp) = interpolant_nodes(itp)


@inline function _shared_headless(itps, pos, a::Int, b::Int)
    xa = _shared_node_source(itps[a])[pos[a]]
    xb = _shared_node_source(itps[b])[pos[b]]

    # Tie-break on interpolant index to give the heap a total ordering.
    return isless(xa, xb) || (!isless(xb, xa) && a < b)
end


@inline function _shared_siftdown!(heap, nheap, root, itps, pos)
    i = root
    k = heap[i]

    @inbounds while true
        left = i << 1
        left > nheap && break

        right = left + 1

        child =
            right <= nheap &&
            _shared_headless(itps, pos, heap[right], heap[left]) ?
            right :
            left

        _shared_headless(itps, pos, heap[child], k) || break

        heap[i] = heap[child]
        i = child
    end

    @inbounds heap[i] = k

    return nothing
end


function _shared_nodes(itps)
    N = length(itps)
    N == 0 && throw(ArgumentError("interpolant collection must be nonempty"))

    first_nodes = _shared_node_source(itps[1])
    X = eltype(first_nodes)

    # pos[k] = current node index for interpolant k
    pos = Vector{Int}(undef, N)

    # Heap contains interpolant indices.
    #
    # The key for stream k is
    #
    #     _shared_node_source(itps[k])[pos[k]]
    #
    heap = Vector{Int}(undef, N)
    nheap = 0

    @inbounds for k in 1:N
        nodes = _shared_node_source(itps[k])

        if !isempty(nodes)
            pos[k] = firstindex(nodes)

            nheap += 1
            heap[nheap] = k
        end
    end

    # No nodes anywhere.
    nheap == 0 && return X[]

    # ------------------------------------------------------------------------
    # Initial heap construction: O(N)
    # ------------------------------------------------------------------------

    for root in (nheap >>> 1):-1:1
        _shared_siftdown!(heap, nheap, root, itps, pos)
    end

    # ------------------------------------------------------------------------
    # Output
    #
    # The union must contain at least as many nodes as the largest individual
    # input, so this sizehint cannot exceed the eventual union length.
    # ------------------------------------------------------------------------

    max_nodes = 0

    @inbounds for k in 1:N
        max_nodes = max(max_nodes, length(_shared_node_source(itps[k])))
    end

    xs = X[]
    sizehint!(xs, max_nodes)

    # Give prev a concrete X type without requiring zero(X).
    @inbounds begin
        k = heap[1]
        prev = _shared_node_source(itps[k])[pos[k]]
    end

    have_prev = false

    # ------------------------------------------------------------------------
    # K-way merge
    # ------------------------------------------------------------------------

    while nheap > 0
        @inbounds begin
            # Heap root is the globally smallest current node.
            k = heap[1]
            nodes = _shared_node_source(itps[k])
            j = pos[k]
            x = nodes[j]

            # Since the merged stream is sorted, duplicates are adjacent.
            if !have_prev || !isequal(x, prev)
                push!(xs, x)
                prev = x
                have_prev = true
            end

            # Advance only the stream we just consumed.
            #
            # Also skip repeated equal nodes within the same interpolant.
            j += 1
            while j <= lastindex(nodes) && isequal(nodes[j], x)
                j += 1
            end

            pos[k] = j

            if j > lastindex(nodes)
                # Stream exhausted. Replace root with final heap leaf.
                heap[1] = heap[nheap]
                nheap -= 1
            end

            # Only the root changed, and its key can only increase.
            # Therefore restoring the heap only requires sift-down.
            if nheap > 0
                _shared_siftdown!(heap, nheap, 1, itps, pos)
            end
        end
    end

    return xs
end


# ============================================================================
# Generic interpolant collection
# ============================================================================

"""
    coerce_to_shared_nodes(itp_collection)

Rebuild every interpolant in the collection on the shared, sorted union of all their nodes, so the
collection is defined on one common node vector and/or type.

Each interpolant's nodes must already be sorted. Works for any backend supplying
[`interpolant_nodes`](@ref) and [`rebuild_interpolant`](@ref), and preserves the container type
(`Vector`, `SVector`, `NTuple`).
"""
function coerce_to_shared_nodes(itp_collection)
    xs = _shared_nodes(itp_collection)

    return map(
        itp -> rebuild_interpolant(itp, xs, itp.(xs)),
        itp_collection,
    )
end


# ============================================================================
# Vector of Fast1DLinearInterpolant
# ============================================================================

function coerce_to_shared_nodes(
    itp_collection::AbstractVector{T},
) where {T<:Fast1DLinearInterpolant}

    xs = _shared_nodes(itp_collection)

    return [
        Fast1DLinearInterpolant(
            xs,
            itp.(xs);
            bc = itp.bc,
            drop_collinear = Val(false),
        )
        for itp in itp_collection
    ]
end


# ============================================================================
# SVector of Fast1DLinearInterpolant
# ============================================================================

function coerce_to_shared_nodes(
    itp_collection::StaticArrays.SVector{N,T},
) where {N,T<:Fast1DLinearInterpolant}

    xs = _shared_nodes(itp_collection)

    return StaticArrays.SVector{N}(
        Fast1DLinearInterpolant(
            xs,
            itp.(xs);
            bc = itp.bc,
            drop_collinear = Val(false),
        )
        for itp in itp_collection
    )
end


# ============================================================================
# NTuple of Fast1DLinearInterpolant
# ============================================================================

function coerce_to_shared_nodes(
    itp_collection::NTuple{N,<:Fast1DLinearInterpolant},
) where {N}

    xs = _shared_nodes(itp_collection)

    return ntuple(
        i -> begin
            itp = itp_collection[i]

            Fast1DLinearInterpolant(
                xs,
                itp.(xs);
                bc = itp.bc,
                drop_collinear = Val(false),
            )
        end,
        Val(N),
    )
end