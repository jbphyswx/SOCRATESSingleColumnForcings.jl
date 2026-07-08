#=
Insertion of surface ("ground") values into the atmospheric column: `combine_air_and_ground_data`
appends/splices a ground slice onto a variable along the vertical dimension, handling scalar,
labeled, and unlabeled ground data and singleton-dimension broadcasting.
=#

"""
    combine_air_and_ground_data

var  : the data variable or string variable name for data in the air
varg : the data variable or string variable name for data on the ground
data : is a container from which the data can be accessed as a string in data[var(g)] form

reshape_ground : expand the ground data to the same size as bottom slice of the air data (so you can pass in for example a single scalar or similar reduced dimension varg)
assume_monotonic : assume that the ground value is actually below everything in the array... in reality sometimes we have for example a ground pressure above that of the minimum in this array... (should be faster than using insert_location )
# might not need this anymore cause if insert_location is integer that jut works

# atlas stated "We use hourly pressure level data interpolated onto a horizontal grid of 0.25° × 0.25° and 37 pressure levels from its native 137 hybrid sigma/pressure levels and 30 km horizontal grid." so i guess sometimes this causes slight problems
"""
function combine_air_and_ground_data(
    var,
    varg,
    concat_dim::Union{Int, String};
    data = nothing,
    reshape_ground::Bool = true,
    insert_location::Union{Symbol, AbstractArray, Function} = :end,
)
    # if var is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata = isa(var, String) ? data[var] : var
    vardatag = isa(varg, String) ? data[varg] : varg
    concat_dim = isa(concat_dim, String) ? dim_num(concat_dim, vardata) : concat_dim # string to numbers

    # handle the ground data alignment
    if reshape_ground
        if isa(vardatag, Number) # allow for you to pass a single number and still concat
            sz_vardatag = collect(size(vardata)) # array
            sz_vardatag[concat_dim] = 1
            vardatag = fill(vardatag, sz_vardatag...) # create array full with just this one value
        elseif hasmethod(NC.dimnames, Tuple{typeof(vardatag)}) # labeled array-like data (CFVariable moved to CommonDataModel in NCDatasets v0.14)
            dimnamesg = NC.dimnames(vardatag)
            # in same order as vardata just in case
            vardatag = reshape(vardatag, (size(vardatag)..., 1)) # add trailing singleton for lev
            dimnamesg = [dimnamesg..., "lev"]
            if hasmethod(NC.dimnames, Tuple{typeof(vardata)})
                dimnames = NC.dimnames(vardata)
                vardatag = permutedims(vardatag, [findfirst(x -> x == dim, dimnames) for dim in dimnamesg]) # permute into the right order
            elseif !isnothing(data) && haskey(data, "T")
                dimnames = NC.dimnames(data["T"])
                vardatag = permutedims(vardatag, [findfirst(x -> x == dim, dimnames) for dim in dimnamesg])
            end
        elseif isa(vardatag, AbstractArray) # an unlabeled array, i think we can't guarantee then that the lev axis exists at all or what existing dimensions are so this just assumes everything is correct except the lev dimension being there
            # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
            sz_vardatag = collect(size(vardatag)) # array
            if ndims(vardatag) == (ndims(vardata) - 1)
                ## TODO USE add_dim
                insert!(sz_vardatag, concat_dim, 1) # insert sz 1 at this location
                vardatag = reshape(vardatag, sz_vardatag...) # reshape (just adds the singleton dimension in)
            elseif ndims(vardatag) != ndims(vardata)
                error(
                    "size mismath, var and varg should have the same number of dimensions or one less (for a missing lev dimension)",
                )
            end
        end
    end

    # Broadcast out singleton dimensions to match sizes except for concat dim.
    # We avoid integer-division based repeat factors because they can silently
    # produce zeros for non-divisible shapes.
    sz_vardata = collect(size(vardata)) # array
    sz_vardatag = collect(size(vardatag)) # array
    length(sz_vardata) == length(sz_vardatag) ||
        error("size mismatch, var and varg must have the same number of dimensions after reshape")
    num_repeat = ones(Int, length(sz_vardata))
    num_repeatg = ones(Int, length(sz_vardata))
    for i in eachindex(sz_vardata)
        i == concat_dim && continue
        sz_a = sz_vardata[i]
        sz_g = sz_vardatag[i]
        if sz_a == sz_g
            continue
        elseif sz_a == 1
            num_repeat[i] = sz_g
        elseif sz_g == 1
            num_repeatg[i] = sz_a
        else
            error("size mismatch in non-concat dimension $(i): var has $(sz_a), varg has $(sz_g)")
        end
    end
    vardatag = repeat(vardatag, num_repeatg...)
    vardata = repeat(vardata, num_repeat...)

    if insert_location == :end # here we assume the ground values are below the values we add automatically. (we used to use identity fcn but i think we don't need that)
        vardata = cat(vardata, vardatag; dims = concat_dim) # we concatenate it at the end...
    elseif isa(insert_location, Function) # use function to determine where to insert our ground values, uses search sorted assuming the original array is already sorted (speedup from binary search)
        mapslice_func = function (vect; by = insert_location) # write this way cause can't define func inside conditional unless anonymous?, see https://github.com/JuliaLang/julia/issues/15602 , https://stackoverflow.com/a/65660721
            vardata = vect[1:(end - 1)]
            vardatag = vect[end]
            return insert_sorted(vardata, vardatag; by = insert_location)
        end
        vardata = cat(vardata, vardatag; dims = concat_dim)
        vardata = mapslices(mapslice_func, vardata; dims = (concat_dim,))

    elseif isa(insert_location, AbstractArray) # use provided indices to determine where to input values...
        if isa(vardatag, Number) # single value, just splice in our value
            vardata = cat(
                selectdim(vardata, concat_dim, 1:(insert_location - 1)),
                vardatag,
                selectdim(vardata, concat_dim, insert_location:size(vardata, concat_dim));
                dims = concat_dim,
            ) # insert the slice there...
        elseif isa(vardatag, AbstractArray) # an unlabeled array, i think we can't guarantee then that the lev axis exists at all or what existing dimensions are so this just assumes everything is correct except the lev dimension being there
            # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
            insert_location_arr = insert_location
            if ndims(insert_location_arr) == (ndims(vardata) - 1)
                sz_insert = collect(size(insert_location_arr))
                insert!(sz_insert, concat_dim, 1)
                insert_location_arr = reshape(insert_location_arr, sz_insert...)
            elseif ndims(insert_location_arr) != ndims(vardata)
                error("insert_location must have the same number of dimensions as vardata (or one less without concat_dim)")
            end

            sz_data = collect(size(vardata))
            sz_ground = collect(size(vardatag))
            sz_insert = collect(size(insert_location_arr))

            target_sz = copy(sz_data)
            for i in eachindex(target_sz)
                if i == concat_dim
                    sz_insert[i] == 1 ||
                        error("insert_location must have singleton concat dimension before insertion")
                    continue
                end
                target_i = max(sz_data[i], sz_ground[i], sz_insert[i])
                ((sz_data[i] == 1) || (sz_data[i] == target_i)) ||
                    error("var shape mismatch in dimension $(i): expected $(target_i) or 1, got $(sz_data[i])")
                ((sz_ground[i] == 1) || (sz_ground[i] == target_i)) ||
                    error("varg shape mismatch in dimension $(i): expected $(target_i) or 1, got $(sz_ground[i])")
                ((sz_insert[i] == 1) || (sz_insert[i] == target_i)) ||
                    error("insert_location shape mismatch in dimension $(i): expected $(target_i) or 1, got $(sz_insert[i])")
                target_sz[i] = target_i
            end

            repeat_data = ones(Int, length(target_sz))
            repeat_ground = ones(Int, length(target_sz))
            repeat_insert = ones(Int, length(target_sz))
            for i in eachindex(target_sz)
                i == concat_dim && continue
                repeat_data[i] = target_sz[i] == sz_data[i] ? 1 : target_sz[i]
                repeat_ground[i] = target_sz[i] == sz_ground[i] ? 1 : target_sz[i]
                repeat_insert[i] = target_sz[i] == sz_insert[i] ? 1 : target_sz[i]
            end

            vardata = repeat(vardata, repeat_data...)
            vardatag = repeat(vardatag, repeat_ground...)
            insert_location_arr = repeat(insert_location_arr, repeat_insert...)

            mapslice_func = function (vect) # write this way cause can't define func inside conditional unless anonymous?, see https://github.com/JuliaLang/julia/issues/15602 , https://stackoverflow.com/a/65660721
                T = promote_type(typeof(vect[1]), typeof(vect[end - 1]))
                vardata = T.(collect(vect[1:(end - 2)]))
                vardatag = T(vect[end - 1])
                insert_location = Int(vect[end]) # undo if was coerced to FT
                # insert (a little more complicated now cause we aren't exactly getting the same size output... dont think that actually matters for mapslices...)
                return insert!(vardata, insert_location, vardatag)
            end
            vardata = cat(vardata, vardatag, insert_location_arr; dims = concat_dim)
            vardata = mapslices(mapslice_func, vardata; dims = (concat_dim,))
        end
    else
        error("unsupported input type for variable insert_location") # catch what would otherwise silently fail and return vardata
    end

    return vardata
end
