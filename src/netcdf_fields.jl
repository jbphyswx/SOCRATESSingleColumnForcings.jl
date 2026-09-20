#=
NetCDF field boundary: materialize lazy NCDatasets arrays into plain Julia `Array`s (replacing
`Missing` with `NaN`) and read profiles/slabs while asserting the expected non-singleton dimensions.
This is the guard that keeps `Missing` and lazy DiskArrays out of the numeric routines downstream.
=#

"""
    read_profile_at_time(vardata, z_dim_num, time_dim_num, time_index)

Read one time slice as a vertical profile. Dimensions other than `z_dim_num`
must be singleton after selecting `time_index`.
"""
function read_profile_at_time(vardata, z_dim_num::Int, time_dim_num::Int, time_index::Int, ::Type{FT} = Base.nonmissingtype(eltype(vardata))) where {FT}
    FTR = (FT <: Nothing) ? Base.nonmissingtype(eltype(vardata)) : FT
    vardata = Array{FTR}(vardata)
    profile = selectdim(vardata, time_dim_num, time_index)
    z_dim_num = z_dim_num - (time_dim_num < z_dim_num)
    _assert_shape_contract(profile, (z_dim_num,), "read_profile_at_time")
    return vec(profile)
end

"""
    read_profiles_over_time(vardata, z_dim_num, time_dim_num; time_indices = :)

Read a variable as a `(z, time)` matrix while asserting that every other
dimension is singleton.
"""
function read_profiles_over_time(
    vardata,
    z_dim_num::Int,
    time_dim_num::Int,
    ::Type{FT} = Base.nonmissingtype(eltype(vardata))
    ;
    time_indices = Colon(),
) where {FT}
    FTR = (FT <: Nothing) ? Base.nonmissingtype(eltype(vardata)) : FT
    vardata = Array{FTR}(vardata)
    slab = time_indices isa Colon ? vardata : selectdim(vardata, time_dim_num, time_indices)
    _assert_shape_contract(slab, (z_dim_num, time_dim_num), "read_profiles_over_time")

    perm = (
        z_dim_num,
        time_dim_num,
        (dim for dim in 1:ndims(slab) if dim != z_dim_num && dim != time_dim_num)...,
    )
    slab = permutedims(slab, perm)
    return reshape(slab, size(slab, 1), size(slab, 2))
end

"""
    _materialize(x)

Materialize any lazy `AbstractArray` (e.g. `DiskArrays.BroadcastDiskArray` produced
by NCDatasets arithmetic into a plain Julia `Array`.
If the element type includes `Missing`, those values are replaced with `NaN` via
`NCDatasets.nomissing(arr, NaN)`, narrowing the element type to `Float64`.

This intentionally mirrors NCDatasets docs usage:
- array conversion: `Array(var)`
- missing handling: `NCDatasets.nomissing(...)`

This is the boundary guard that prevents `Missing` from ever propagating into the
package's interpolation and computation routines.
"""
function _materialize(x::AbstractArray{T}, ::Type{FT} = Base.nonmissingtype(eltype(x))) where {T, FT}
    FTR = (FT <: Nothing) ? Base.nonmissingtype(eltype(x)) : FT
    arr = (x isa Array) ? x : Array{T}(x)
    if Missing <: eltype(x)
        FNMT = Base.nonmissingtype(eltype(arr))
        out = NC.nomissing(arr, FNMT(NaN))
        out = (FNMT === FTR) ? out : convert(Array{FTR}, out)
    else
        out = (FTR === T) ? arr : convert(Array{FTR}, arr)
    end
    return out
end

"""
    _materialize!(out, x)

In-place [`_materialize`](@ref): load `x` into the caller's `out` and return the materialized array.
When `T` includes `Missing` the `Missing`-stripped result is a fresh array (`nomissing` narrows the
element type), so use the returned value rather than assuming `out` holds it.
"""
function _materialize!(out::AbstractArray{T}, x::AbstractArray{FT}) where {T, FT}
    NCDatasets.load!(x, out, fill(Colon(), ndims(x))...)
    if Missing <: T
        FNMT = Base.nonmissingtype(T)
        out = NC.nomissing(out, FNMT(NaN))
        out = (FNMT === FTR) ? out : convert(Array{FTR}, out)
    end
    return out
end
