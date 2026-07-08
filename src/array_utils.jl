#=
Array / dimension utilities: axis manipulation, dimension lookup by name, sorted insertion, and a
small NamedTuple accessor. No dependency on the interpolation, thermodynamics, or NetCDF layers
beyond `NCDatasets.dimnames` for label-based dimension lookup.
=#

"""
    swapaxes(a, dim1, dim2)

Swaps the two dimensions of an array.
"""
function swapaxes(a, dim1, dim2)
    perm = collect(1:ndims(a))
    perm[dim1] = dim2
    perm[dim2] = dim1
    return permutedims(a, perm)
end

"""
    align_along_dimension(v,dim)

Assumes v has one non-singleton dimension which will be aligned along dim.
if dim does not exist, v is expanded till dim exists
"""
function align_along_dimension(v, dim)

    sz = size(v)
    nd = ndims(v)
    if nd < dim # add dims if they don't exist
        v = reshape(v, sz..., fill(1, dim - nd)...)
    end

    # get nonsingleton dimension
    nonsingleton_ind = findall(sz .> 1)
    l_ns = length(nonsingleton_ind)
    if l_ns > 1
        error("v has more than one non-singleton dimension")
    elseif l_ns < 1
        error("v has no non-singleton dimension")
    end

    v = swapaxes(v, nonsingleton_ind[1], dim) # move the non-singleton dimension to the dim position
    return v
end

"""
    add_dim(vardata, dimnum)

Allows you to add a dimension to a variable at the specified location
"""
function add_dim(vardata, dimnum)
    # both data need to be labeled
    sz_vardata = collect(size(vardata)) # array
    insert!(sz_vardata, dimnum, 1) # insert sz 1 at this location
    return reshape(vardata, sz_vardata...) # reshape
end

function _assert_shape_contract(vardata, keep_dims::Tuple, name::AbstractString)
    unexpected_dims = [
        dim for dim in 1:ndims(vardata) if size(vardata, dim) > 1 && !(dim in keep_dims)
    ]
    isempty(unexpected_dims) || error(
        "$name expected only dims $(collect(keep_dims)) to be non-singleton, got size $(size(vardata)) with extra non-singleton dims $(unexpected_dims)",
    )
    return nothing
end

function insert_sorted(vect, val; by = identity, rev = false)
    index = searchsortedfirst(vect, val; by = by, rev = rev) #find index at which to insert x
    return insert!(vect, index, val) #insert x at index
end

"""
    dim_num

get the dimension number of dim from nc_data
"""
dim_num(dim::Number, nc_data) = dim
function dim_num(dim::String, nc_data)
    if hasmethod(NC.dimnames, Tuple{typeof(nc_data)})
        dimnames = NC.dimnames(nc_data)
        dim_num = findfirst(x -> x == dim, dimnames)
        return dim_num
    elseif isa(nc_data, NC.NCDataset)
        error("dimension number for a dataset is not well defined, use a specific variable with dimension labels instead")
    elseif isa(nc_data, AbstractArray)
        error("cannot find dimension $(dim) in unlabeled data nc_data, pass in a labeled NCDatasets variable instead...")
    else
        error("unsupported input type for nc_data $(typeof(nc_data))")
    end
end

set_property(x::NamedTuple, field::Symbol, value) = merge(x, (field => value,))  # edit the value in a named tuple
