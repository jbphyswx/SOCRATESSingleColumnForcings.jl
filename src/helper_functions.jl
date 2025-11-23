#= ...
created jan 4 2023
@jbphyswx
... =#

##
# - Maybe swith to using DimensionalData or something so we can have derived quantities (i.e. ts <thermodynamic state>) as labeled arrays and not have to keep passing their parents around...
# -- seems to be a little cumbersome and no easy way to do netCDF reads and conversion to array/dataset like types
##

# ================================================================================================================================================================================ #
include("interpolating_methods.jl") # 
include("les_reader_helper.jl") # implements get_data_new_z_t_LES() for LES data instead of the input forcing files from Atlas
# ================================================================================================================================================================================ #

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

"""
    data_to_tsg(data; thermo_params)

Does this cause we don't have surface q values as of rn
"""
function data_to_tsg(data; thermo_params::TDPS)

    Tg = data["Tg"]
    pg = data["Ps"]

    # not sure what to do at surface, so assuming saturation at surface
    pvg = TD.saturation_vapor_pressure.(thermo_params, Tg, TD.Liquid())
    molmass_ratio = TDP.molmass_ratio(thermo_params)

    qg = (1 / molmass_ratio) .* pvg ./ (pg .- pvg) #Total water mixing ratio at surface , assuming saturation

    tsg = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qg)
    return tsg
end

"""
    data_to_ts
"""
function data_to_ts(data; do_combine_air_and_ground_data::Bool = false, thermo_params::TDPS) # add type data is ncdataset

    T = data["T"]
    q = data["q"]
    p = data["lev"]

    # put p/lev on correct dimension...
    p = align_along_dimension(p, get_dim_num("lev", data["T"]))

    ts = TD.PhaseEquil_pTq.(thermo_params, p, T, q)

    if do_combine_air_and_ground_data
        concat_dim = get_dim_num("lev", T) # phaseequil reduces us down to lev... can't seem to just apply ufunc...
        tsg = data_to_tsg(data; thermo_params)
        ts = combine_air_and_ground_data(ts, tsg, concat_dim)
    end

    return ts
end

function insert_sorted(vect, val; by = identity, rev = false)
    index = searchsortedfirst(vect, val; by = by, rev = rev) #find index at which to insert x
    return insert!(vect, index, val) #insert x at index
end







function lev_to_z_from_LES_output_column(
    tsz::AbstractArray{TS},
    lesz,
    lesp;
    thermo_params::TDPS,
    interp_method::AbstractInterpolationMethod = FastLinear1DInterpolation,
) where {TS <: TDTS}
    # Load pressure from the thermodynamic state
    p_in = TD.air_pressure.(thermo_params, tsz)
    # Interpolate z from the LES output to the pressure levels of the thermodynamic state
    z_out = interpolate_1d(p_in, lesp, lesz, interp_method; bc = ErrorBoundaryCondition())
    return z_out
end

"""
We're trying to use the LES output p, z to get the input p's z's since using the thickness equation may have meant we're off a little
However, interpolating back to the input data's scale could be bad and lose the resolution gains we have... we'll have to see how it goes...
We only have the forcing data at that scale so...
 
"""
function lev_to_z_from_LES_output(
    ts,
    tsg;
    thermo_params::TDPS,
    data,
    assume_monotonic::Bool = false,
    flight_number::Int,
    forcing_type::Symbol,
    ground_indices = nothing,
)

    dimnames = NC.dimnames(data["T"]) # use this as default cause calculating ts doesn't maintain dim labellings
    lev_dim_num = findfirst(x -> x == "lev", dimnames)
    ldn = lev_dim_num
    L = size(ts, lev_dim_num)

    time_dim_num = findfirst(x -> x == "time", dimnames)
    tdn = time_dim_num
    Lt = size(ts, time_dim_num)

    if !assume_monotonic
        LES_data = open_atlas_les_output(flight_number, forcing_type)[forcing_type]
        p_LES = LES_data["PRES"] # 2D,time varying
        z_LES = LES_data["z"] # 1D, constant
        # LES_min_z, LES_max_z = extrema(z_LES)

        dimnames_LES = NC.dimnames(p_LES) # use this as default cause calculating ts doesn't maintain dim labellings (should be time, z)
        # ldn_LES =  findfirst(x -> x == "z", dimnames_LES)
        tdn_LES = findfirst(x -> x == "time", dimnames_LES)

        # handle times...

        # overall
        summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
        SOCRATES_summary = NC.Dataset(summary_file, "r")
        flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
        initial_time = SOCRATES_summary["reference_time"][flight_ind] - Dates.Hour(12) # simulation start time is 12 hours before reference time, but set up socrates summary to give them all the same reference... (is in jupyter notebook somewhere)

        # input ( i think this uses same calendar as overall? idk...)
        t_in = data["tsec"] # should just be seconds past initial_time
        t_base = Dates.DateTime(string(data["bdate"][:]), Dates.DateFormat("yymmdd")) + Dates.Year(2000) # the base Date (using bdate not nbdate cause nbdate seems to have  bug in flight 9 (extra 0 in month spot))
        t_in = t_base .+ Dates.Second.(t_in) # the actual dates

        t_in = (t_in .- initial_time) ./ Dates.Millisecond(1) # make it relative to the start of the simulation
        t_in ./= 1000 # Milliseconds to seconds

        # LES (this has no calendar, is just t in seconds, i *think* it's on the same t as 
        t_les = LES_data["time"] # should just be seconds past initial_time
        t_les = (t_les .- t_les[1]) * (24 * 3600 * 1000) # make it relative to the start of the simulation, and convert days to miliseconds (then round, bc Dates can't handle non integer amounts)
        t_les ./= 1000 # Milliseconds to seconds
        t_les = NCDatasets.nomissing(t_les[:])


        p_LES = NC.nomissing(p_LES[:]) .* 100 # convert to hPa to Pa
        z_LES = NC.nomissing(z_LES[:])

        # interpolate les output data to input times (instead of just choosing the `closest` hour like we did before... (need to apply by row)
        # i think you need to map this bc of the splines... but there's probably some way...
        p_LES = mapslices(
            x -> interpolate_1d(t_in, t_les, x, FastLinear1DInterpolation; bc = NearestBoundaryCondition()),
            p_LES;
            dims = 2,
        ) # apply to each row, use nearest interp but outside les time bounds is bogus

        tsz = ts

        FTLES = eltype(z_LES)

        # preallocate z
        s_tsz = collect(size(ts))
        s_tsz[ldn] += 1 # we just want to add a single slice of zeros for the ground -- these are all ocean cases so this should be fine for now...
        z = Array{FTLES}(undef, s_tsz...) # should be same size as ts

        # iterate over columns in z (ldn)

        # One problem is the LES doesn't span the full range of the input... so should we do the LES for the LES range and then use the thickness equation outside that?
        increasing_p = false

        # add one for ground
        len = length(selectdim(p_LES, tdn_LES, 1)[:])  # number of elements in the slice
        p_LES_t = Vector{FTLES}(undef, len)

        len = length(selectdim(tsz, tdn, 1)[:])
        tsz_t_no_g = Vector{eltype(tsz)}(undef, len) # preallocate for type stability
        tsz_t = Vector{eltype(tsz)}(undef, len + 1)       # preallocate for type stability
        p_in_t_no_g = Vector{FTLES}(undef, len) # preallocate for type stability
        p_in_t = Vector{FTLES}(undef, len + 1)       # preallocate for type stability


        for i_t in 1:Lt

            # LES does not span the full range of the input, so we need to do the LES for the LES range and then use the thickness equation outside that
            # Get min/max inds for input data that are covered by LES data

            tsz_t_no_g .= selectdim(tsz, tdn, i_t)[:]
            tsgz = tsg[i_t]

            p_in_t_no_g .= TD.air_pressure.(thermo_params, tsz_t_no_g[:])
            # p_in_min, p_in_max = extrema(p_in_t_no_g)

            p_LES_t .= selectdim(p_LES, tdn_LES, i_t)[:]
            p_LES_min, p_LES_max = extrema(p_LES_t)

            p_s_in = TD.air_pressure(thermo_params, tsgz)

            if isnothing(ground_indices)
                index = searchsortedfirst(tsz_t_no_g, tsgz; by = x -> TD.air_pressure(thermo_params, x), rev = false) # find where the ground value would be inserted...
            else
                index = selectdim(ground_indices, tdn, i_t)[:][] # should just be one item
            end

            # sort them all same direction (LES to match input, and we're going with increasing_p for simplifying logic below (decreasing z))
            increasing_p = p_in_t_no_g[1] < p_in_t_no_g[end] # if the first pressure is lower than the last, we're increasing
            if increasing_p
                sort!(p_LES_t)
                sort!(z_LES, rev = true)
            else
                reverse!(p_in_t_no_g)
                reverse!(tsz_t_no_g)
                sort!(p_LES_t)
                sort!(z_LES, rev = true)
                index = length(p_in_t) - index + 1 # reverse the index
            end

            # insert!(tsz_t, index, tsgz) # only seems to be an in place option...
            # insert!(p_in_t, index, p_s_in) # only seems to be an in place option...
            tsz_t = [tsz_t_no_g[1:(index - 1)]; tsgz; tsz_t_no_g[index:end]]
            p_in_t = [p_in_t_no_g[1:(index - 1)]; p_s_in; p_in_t_no_g[index:end]]

            # (we're going increasing p direction) 
            i_t_min_p = findfirst(x -> x > p_LES_min, p_in_t) # first ind with pressure higher than the lowest pressure in the LES data (going increasing_p direction)
            i_t_max_p = findlast(x -> x < p_LES_max, p_in_t)  #  last ind with pressure lower than the highest pressure in the LES data (going increasing_p direction)

            squeeze(a) = dropdims(a, dims = tuple(findall(size(a) .== 1)...))

            # handle data covered by LES
            selectdim(z, tdn, i_t)[i_t_min_p:i_t_max_p] =
                lev_to_z_from_LES_output_column(tsz_t[i_t_min_p:i_t_max_p], z_LES[:], p_LES_t; thermo_params)


            # Handle z's lower than LES Data (highest pressure to sfc) | use the thickness equation
            new_z = lev_to_z_column(tsz_t[i_t_max_p:end]; thermo_params)
            new_dz = new_z[2:end] .- new_z[1:(end - 1)] # get the dz from the thickness equation, but instead of going up from the ground, we're gonna flip it to go down from the last good z
            new_dz = cumsum(new_dz) # sum up our dz
            new_z = selectdim(z, tdn, i_t)[i_t_max_p] .+ new_dz
            selectdim(z, tdn, i_t)[(i_t_max_p + 1):end] = new_z



            # Handle z's higher than LES Data | use the thickness equation and add to the top of the LES data
            selectdim(z, tdn, i_t)[1:(i_t_min_p - 1)] =
                lev_to_z_column(tsz_t[1:i_t_min_p]; thermo_params)[1:(end - 1)] .+ z[i_t_min_p] # go from start to top and add to existing top...


            selectdim(z, tdn, i_t)[:] .-= selectdim(z, tdn, i_t)[index] # subtract out the ground value
        end

        if !increasing_p
            # reverse!(z, dims = ldn) # put back the way it was
            z = reverse(z, dims = ldn) # put back the way it was
        end

        return z
    else
        error("not implemented")
    end
end


"""
Because we might not have monotonic data, we need to be able to both insert the ground data into our lev array where applicable and calculate our dz accordingly...
I guess in principle you could throw out p < ps but then you wouldn't be able to return an even array out so either way this function would need to exist for padding and such...

tsz should be a one dimensional array consising of [ts..., tg]

# to do -- add capability to use precomputed indices (insert_location)
"""
function lev_to_z_column(tsz; thermo_params::TDPS)

    ts = tsz[1:(end - 1)] # need to make it a vector I guess... (not sure if this screws wit output shape)
    tsg = tsz[end]
    R_d = TDP.R_d(thermo_params)
    grav = TDP.grav(thermo_params)

    L = length(ts)

    index = searchsortedfirst(ts, tsg; by = x -> TD.air_pressure(thermo_params, x), rev = false) # find where the ground value would be inserted...
    insert!(ts, index, tsg) # only seems to be an in place option...
    tsz = ts # replace with the reordered version

    Tvz = TD.virtual_temperature.(thermo_params, tsz) # virtual temp, we havent returned these for now...
    pz = TD.air_pressure.(thermo_params, tsz)
    Lz = L + 1 # cause we extended it using the ground...
    Tv_bar = Statistics.mean((Tvz[1:(Lz - 1)], Tvz[2:Lz]))
    p_frac = pz[2:Lz] ./ pz[1:(Lz - 1)]
    dz = @. (R_d * Tv_bar / grav) * log(p_frac)

    # sum up from the bottom then subtract the height of the ground
    # z = reverse(cumsum(reverse(dz))) #cumsum(dz) # grid is already defined from  (from Grid.jl)
    # z = [z..., 0] # is this the right order? seems so based on the cat below but idk... if so might have to flip index to be L - index or something like that? (seems to be so...)

    # Preallocate output
    L = length(dz)
    # z = zeros(eltype(dz), L+1)  # last element = ground
    z = similar(dz, L + 1)  # last element = ground
    z[end] = eltype(dz)(0)  # set ground level to 0
    for i in L:-1:1 # # Fill from top down
        z[i] = z[i + 1] + dz[i]
    end

    # after the padding, z for the ground should be at the right index, and we can just subtract it out from the array...
    z = z .- z[index]
    return z
end

# convert pressure to altitude...
"""
    lev_to_z

TODO: document
"""
function lev_to_z(p::FT, T::FT, q::FT, pg::FT, Tg::FT, qg::FT; thermo_params::TDPS, data) where {FT <: Real}
    ts = TD.PhaseEquil_pTq.(thermo_params, p, T, q)
    tsg = TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qg)
    return lev_to_z(ts, tsg; data = data, thermo_params)
end


"""
    lev_to_z

ts is thermodynamic state
tsg is thermodynamic state for ground

if assume monotonic, everything should already be in the right order and we can use the vectorized version, otherwise we will use lev_to_z column applied column by column w/ mapslices
"""
function lev_to_z(
    ts::AbstractArray{TS},
    tsg::AbstractArray{TS};
    thermo_params::TDPS,
    data,
    assume_monotonic::Bool = false,
) where {TS <: TDTS}

    dimnames = NC.dimnames(data["T"]) # use this as default cause calculating ts doesn't maintain dim labellings
    lev_dim_num = findfirst(x -> x == "lev", dimnames)
    ldn = lev_dim_num
    L = size(ts, lev_dim_num)

    if !assume_monotonic
        tsz = combine_air_and_ground_data(ts, tsg, ldn; data, reshape_ground = true, insert_location = :end) # here we just want to assume monotonic so we can pass those to this fcn
        z = mapslices(x -> lev_to_z_column(x; thermo_params), tsz; dims = ldn) # need to make a stack cause that's all mapslices can take...
    else

        R_d = TDP.R_d(thermo_params) # TD.Parameters.R_d(thermo_params)
        grav = TDP.grav(thermo_params)  # TD.Parameters.grav(thermo_params)


        tsz = combine_air_and_ground_data(
            ts,
            tsg,
            ldn;
            data = data,
            reshape_ground = true,
            insert_location = x -> TD.air_pressure(thermo_params, x),
        ) # this doesnt work in reality cause sometimes ps is more than the lowest value in lev... so we use the not assume_monotonic version mostly

        #actually we need to split dz into pos or neg depending on whether or not it's above ground... maybe best just to have fcn that goes along each column...

        Tvz = TD.virtual_temperature.(thermo_params, tsz) # virtual temp, we havent returned these for now...
        pz = TD.air_pressure.(thermo_params, tsz)
        Lz = L + 1 # cause we extended it using the ground...
        Tv_bar = Statistics.mean((selectdim(Tvz, ldn, 1:(Lz - 1)), selectdim(Tvz, lev_dim_num, 2:Lz)))
        p_frac = selectdim(pz, lev_dim_num, 2:Lz) ./ selectdim(pz, ldn, 1:(Lz - 1))
        dz = @. (R_d * Tv_bar / grav) * log(p_frac)

        z = reverse(cumsum(reverse(dz; dims = ldn); dims = ldn); dims = ldn) #cumsum(dz) # grid is already defined from  (from Grid.jl)
        s_sz = collect(size(z)) # should now be same dims as T
        s_sz[ldn] = 1 # we just want to add a single slice of zeros for the ground -- these are all ocean cases so this should be fine for now...
        _z0 = zeros(Float64, s_sz...)
        z = cat(z, _z0; dims = ldn)
    end
    return z
end

function z_from_data(data; thermo_params::TDPS)
    ts = data_to_ts(data; thermo_params, do_combine_air_and_ground_data = false)
    tsg = data_to_tsg(data; thermo_params)
    return lev_to_z(ts, tsg; thermo_params, data)
end


"""
    get_ground_insertion_indices

Get the indices where the ground tsg would fit into the array ts...
"""
function get_ground_insertion_indices(
    ts::AbstractArray{TS},
    tsg::AbstractArray{TS},
    concat_dim::Union{Int, String};
    thermo_params::TDPS,
    data,
) where {TS <: TDTS}
    function mapslice_func(vect; thermo_params = thermo_params, by = x -> TD.air_pressure(thermo_params, x))
        vardata = vect[1:(end - 1)]
        vardatag = vect[end]
        index = searchsortedfirst(vardata, vardatag; by = by, rev = false)
        return index
    end
    # reshape ( TODO!! # use add_dim )
    sz_tsg = collect(size(tsg)) # array
    insert!(sz_tsg, concat_dim, 1) # insert sz 1 at this location
    tsg = reshape(tsg, sz_tsg...) # reshape (just adds the singleton dimension in)
    # concat and calculate
    ts = cat(ts, tsg; dims = concat_dim)
    return mapslices(mapslice_func, ts; dims = [concat_dim])
end




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
    concat_dim = isa(concat_dim, String) ? get_dim_num(concat_dim, vardata) : concat_dim # string to numbers

    # handle the ground data alignment
    if reshape_ground
        if isa(vardatag, Number) # allow for you to pass a single number and still concat
            sz_vardatag = collect(size(vardata)) # array
            sz_vardatag[concat_dim] = 1
            vardatag = fill(vardatag, sz_vardatag...) # create array full with just this one value
        elseif isa(vardatag, NC.CFVariable) # vardatag should not have lev as a dimension so we need to add it in (and I guess check the others are in the same order)
            dimnamesg = NC.dimnames(vardatag)
            dimnames = NC.dimnames(vardata)
            # in same order as vardata just in case
            vardatag = reshape(vardatag, (size(vardatag)..., 1)) # add trailing singleton for lev
            dimnamesg = [dimnamesg..., "lev"]
            vardatag = permutedims(vardatag, [findfirst(x -> x == dim, dimnames) for dim in dimnamesg]) # permute into the right order
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

    # broadcast out to match sizes except for concat dim
    sz_vardata = collect(size(vardata)) # array
    sz_vardatag = collect(size(vardatag)) # array
    num_repeat = sz_vardatag .÷ sz_vardata # integer division
    num_repeatg = sz_vardata .÷ sz_vardata # integer division
    num_repeat[concat_dim] = 1 # don't repeat along concat dim
    num_repeatg[concat_dim] = 1 # don't repeat along concat dim
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
            mapslice_func = function (vect) # write this way cause can't define func inside conditional unless anonymous?, see https://github.com/JuliaLang/julia/issues/15602 , https://stackoverflow.com/a/65660721
                vardata = vect[1:(end - 2)]
                vardatag = vect[end - 1]
                insert_location = Int(vect[end]) # undo if was coerced to FT
                # insert (a little more complicated now cause we aren't exactly getting the same size output... dont think that actually matters for mapslices...)
                return insert!(vardata, insert_location, vardatag)
            end
            vardata = cat(vardata, vardatag, insert_location; dims = concat_dim)
            vardata = mapslices(mapslice_func, vardata; dims = (concat_dim,))
        end
    else
        error("unsupported input type for variable insert_location") # catch what would otherwise silently fail and return vardata
    end

    return vardata
end




"""
    get_dim_num

get the dimension number of dim from nc_data
"""
get_dim_num(dim::Number, nc_data) = dim
function get_dim_num(dim::String, nc_data::NC.CFVariable)
    dimnames = NC.dimnames(nc_data)
    dim_num = findfirst(x -> x == dim, dimnames)
    return dim_num
end
get_dim_num(dim::String, nc_data::NC.NCDataset) =
    error("dimension number for a dataset is not well defined, use a speciifc NC.CFVariable instead")
get_dim_num(dim::String, nc_data::AbstractArray) =
    error("cannot find dimension $(dim) in unlabeled data nc_data, pass in a labeled NCDataset instead...")
get_dim_num(dim::String, nc_data) = error("unsupported input type for nc_data $(typeof(nc_data))")

# function get_dim_num(dim, nc_data = nothing)

#     # convert interp_dim to a number stored in interp_dim_num
#     if isa(dim, Number) # already a string, just returns...
#         dim_num = dim
#     elseif isa(dim, String)
#         if isa(nc_data, NC.CFVariable) # allow finding the number from an ncdataset and a string
#             dimnames = NC.dimnames(nc_data)
#             dim_num = findfirst(x -> x == dim, dimnames)
#         elseif isa(nc_data, NC.NCDataset)
#             @warn("dimension number for a dataset is not well defined, use a speciifc NC.CFVariable instead")
#         elseif isa(nc_data, AbstractArray)
#             error("cannot find dimension $(dim) in unlabeled data nc_data, pass in a labeled NCDataset instead...")
#         else
#             error("unsupported input type for nc_data $(typeof(nc_data))")
#         end
#     end
#     return dim_num
# end



"""
    interp_along_dim

interpolation data

# var : the data to be interpolated
# interp_dim: the name or dimension number along which we will do interpolation

# interp_dim_in: the coordinate on which the data is currently aligned
# interp_dim_out: the coordinate on which we would like to interpolate the data to

the function will create splines along interp_dim_in...
- if we set interp_dim_out, we actually evaluate the function along interp_dim at locations interp_dim_out,otherwise, we just return our unevaluated spline functions

- data_func is a func to be applied to the raw data before it is processed, though perhaps it it most useful if interp_dim_out is unset and we are returning functions...
- vectorize_in means your input is an array and you need to loop over it too (as opposed to just being a fixed template vector)

Note the conservative_regirdder() returns actual output arrays -- it won't work if you are just creating a spline fcn

To Do : decide types
"""
function interp_along_dim(
    var,
    interp_dim::Union{String, Int},
    interp_dim_in;
    interp_dim_out = nothing,
    data = nothing,
    data_func::Union{Function, Nothing} = nothing,
    interp_dim_in_is_full_array::Bool = true,
    # reshape_ground::Bool = true,
    interp_method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    interp_kwargs::Dict = Dict{Symbol, Any}(),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::DCIKT = default_conservative_interp_kwargs,
    A::Union{Nothing, AbstractArray} = nothing,
    Af::Union{Nothing, AbstractArray} = nothing,
    use_svectors::Bool = false,
)
    # would use kwargs but doesn't play nice w/ ODE solver for some reason... instead we get out the parameters we want explicitly and pass them all the time (splatting gives a union typle type object that apparently can't be handled?)
    f_enhancement_factor = get(interp_kwargs, :f_enhancement_factor, 1) # default to 1.0
    f_p_enhancement_factor = get(interp_kwargs, :f_p_enhancement_factor, 1) # default to 1.0
    if conservative_interp
        bc::AbstractBoundaryCondition = create_bc(get(interp_kwargs, :bc, ExtrapolateBoundaryCondition())) # default to extrapolate
    else
        bc = create_bc(get(interp_kwargs, :bc, ErrorBoundaryCondition())) # default to error
    end
    k = get(interp_kwargs, :k, 1) # default to 3



    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata = isa(var, String) ? data[var] : var


    # get the interpolation dimension and combine air and ground data
    interp_dim_num = get_dim_num(interp_dim, vardata)
    # vardata        = combine_air_and_ground_data(vardata ,vardatag, interp_dim_num; data=nothing, reshape_ground=true) # combine air and ground data... (also resolves strings to data)

    if !isnothing(data_func) # apply data_func if we need to
        vardata = data_func(vardata)
    end

    # :: interp_func(x, xp, fp;...) :: #
    if conservative_interp
        interp_func =
            (args...; kwargs...) ->
                conservative_regridder(args...; A = A, Af = Af, conservative_interp_kwargs..., kwargs...)
    else
        interp_func = interpolate_1d
    end
    conversion_func = use_svectors ? create_svector : identity

    # mapslices to apply along timedim, see https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices
    if !interp_dim_in_is_full_array
        if isnothing(interp_dim_out)
            out = mapslices(
                let interp_dim_in = interp_dim_in,
                    interp_method = interp_method,
                    f_enhancement_factor = f_enhancement_factor,
                    f_p_enhancement_factor = f_p_enhancement_factor,
                    bc = bc,
                    k = k
                    # let block for performance of captured variables

                    d -> let d = d
                        # dd -> interp_func(
                        #     dd,
                        #     conversion_func(interp_dim_in),
                        #     conversion_func(d);
                        #     method = interp_method,
                        #     f_enhancement_factor = f_enhancement_factor,
                        #     f_p_enhancement_factor = f_p_enhancement_factor,
                        #     bc = bc,
                        #     k = k,
                        # )
                        build_spline(
                            interp_method,
                            conversion_func(interp_dim_in),
                            conversion_func(d);
                            f_enhancement_factor = f_enhancement_factor,
                            f_p_enhancement_factor = f_p_enhancement_factor,
                            bc = bc,
                            k = k,
                            drop_collinear = false, # messes w/ mapslices inference
                        ) # is already a fcn
                    end
                end,
                vardata,
                dims = (interp_dim_num,),
            ) # will return a lambda fcn that can be evaluated along that dimension. It's here that the svector could be meaningful

            if (drop_collinear = true) &&
               (conversion_func == create_svector) &&
               (interp_method <: FastLinear1DInterpolationMethod) # try to do the reduction in post
                out_type = (SSCF.Fast1DLinearInterpolant{<:StaticArrays.SVector{N, Float64}} where {N}) # broader
                # create empty array of correct type
                out_new = Array{out_type}(undef, size(out)...)
                for i in eachindex(out)
                    out_new[i] = Fast1DLinearInterpolant(out[i].xp, out[i].fp; bc = out[i].bc, drop_collinear = true)
                end
                return out_new
            end

            return out

        else
            out = mapslices(
                let interp_dim_out = interp_dim_out,
                    interp_dim_in = interp_dim_in,
                    interp_method = interp_method,
                    f_enhancement_factor = f_enhancement_factor,
                    f_p_enhancement_factor = f_p_enhancement_factor,
                    bc = bc,
                    k = k
                    # let block for performance of captured variables
                    d -> interp_func(
                        interp_dim_out,
                        conversion_func(interp_dim_in),
                        conversion_func(d);
                        method = interp_method,
                        f_enhancement_factor = f_enhancement_factor,
                        f_p_enhancement_factor = f_p_enhancement_factor,
                        bc = bc,
                        k = k,
                    )
                end,
                vardata,
                dims = (interp_dim_num,),
            ) # lambda fcn will evaluate, svector less meaningful

            if (drop_collinear = true) &&
               (conversion_func == create_svector) &&
               (interp_method <: FastLinear1DInterpolationMethod) # try to do the reduction in post
                out_type = (SSCF.Fast1DLinearInterpolant{<:StaticArrays.SVector{N, Float64}} where {N}) # broader
                # create empty array of correct type
                out_new = Array{out_type}(undef, size(out)...)
                for i in eachindex(out)
                    out_new[i] = Fast1DLinearInterpolant(out[i].xp, out[i].fp; bc = out[i].bc, drop_collinear = true)
                end
                return out_new
            end

            return out
        end
    else # vectorize over input dim values as well as data (no support for vectorize over output dim yet)
        # stack on new catd dimension, then split apart inside the fcn call
        catd = ndims(vardata) + 1
        _input = cat(interp_dim_in, vardata; dims = catd) # although maybe the input is just a vector in which case this won't work..... we could just pass it in
        if isnothing(interp_dim_out)
            out = dropdims(
                mapslices(
                    let interp_method = interp_method,
                        f_enhancement_factor = f_enhancement_factor,
                        f_p_enhancement_factor = f_p_enhancement_factor,
                        bc = bc,
                        k = k
                        # let block for performance of captured variables
                        d -> let d = d
                            # dd -> interp_func(
                            #     dd,
                            #     conversion_func(d[:, 1]),
                            #     conversion_func(d[:, 2]);
                            #     method = interp_method,
                            #     f_enhancement_factor = f_enhancement_factor,
                            #     f_p_enhancement_factor = f_p_enhancement_factor,
                            #     bc = bc,
                            #     k = k,
                            # )
                            build_spline(
                                interp_method,
                                conversion_func(d[:, 1]),
                                conversion_func(d[:, 2]);
                                f_enhancement_factor = f_enhancement_factor,
                                f_p_enhancement_factor = f_p_enhancement_factor,
                                bc = bc,
                                k = k,
                            ) # is already a fcn
                        end
                    end,
                    _input,
                    dims = (interp_dim_num, catd),
                );
                dims = catd,
            ) # will return a lambda fcn that can be evaluated along that dimension, here the svector could be meaningful

            if (drop_collinear = true) &&
               (conversion_func == create_svector) &&
               (interp_method <: FastLinear1DInterpolationMethod) # try to do the reduction in post
                out_type = (SSCF.Fast1DLinearInterpolant{<:StaticArrays.SVector{N, Float64}} where {N}) # broader
                # create empty array of correct type
                out_new = Array{out_type}(undef, size(out)...)
                for i in eachindex(out)
                    out_new[i] = Fast1DLinearInterpolant(out[i].xp, out[i].fp; bc = out[i].bc, drop_collinear = true)
                end
                return out_new
            end

            return out
        else
            out = dropdims(
                mapslices(
                    let interp_dim_out = interp_dim_out,
                        interp_method = interp_method,
                        f_enhancement_factor = f_enhancement_factor,
                        f_p_enhancement_factor = f_p_enhancement_factor,
                        bc = bc,
                        k = k
                        # let block for performance of captured variables
                        d -> interp_func(
                            interp_dim_out,
                            conversion_func(d[:, 1]),
                            conversion_func(d[:, 2]);
                            method = interp_method,
                            f_enhancement_factor = f_enhancement_factor,
                            f_p_enhancement_factor = f_p_enhancement_factor,
                            bc = bc,
                            k = k,
                        )
                    end,
                    _input,
                    dims = (interp_dim_num, catd),
                );
                dims = catd,
            ) # lambda fcn will evaluate

            if (drop_collinear = true) &&
               (conversion_func == create_svector) &&
               (interp_method <: FastLinear1DInterpolationMethod) # try to do the reduction in post
                out_type = (SSCF.Fast1DLinearInterpolant{<:StaticArrays.SVector{N, Float64}} where {N}) # broader
                # create empty array of correct type
                out_new = Array{out_type}(undef, size(out)...)
                for i in eachindex(out)
                    out_new[i] = Fast1DLinearInterpolant(out[i].xp, out[i].fp; bc = out[i].bc, drop_collinear = true)
                end
                return out_new
            end

            return out
        end
    end
end


##########################
#### -- note is this redundant? i think techincally we're create and instantateously evaluate the spline over and over -- maybe we want to just store the spline and then evaluate it at different points?
##########################



# function var_to_new_z(;var=var, varg=nothing, old_z=z, new_z=default_new_z, interp_dim=ldn, data_func=x->reverse(x;dims=ldn), kwargs...)
#     """
#     data interpolation to new z coordinate
#     """

#     interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension

#     # make sure our data is numeric and combine air and ground if necessary...
#     vardata    = isa(var, String) ? data[var ] : var
#     if !isnothing(varg)
#         vardatag   = isa(varg,String) ? data[varg] : varg
#         vardata    = combine_air_and_ground_data(vardata ,vardatag, ldn; reshape_ground=true) # combine air and ground data... (also resolves strings to data)
#     end

#     return interp_along_dim(vardata, interp_dim, reverse(z; dims=interp_dim); interp_dim_out=new_z, data_func=data_func, kwargs...)
# end

"""
    var_to_new_coord

if coord_new is nothing, then will return functions...
"""
function var_to_new_coord(
    var,
    coord_in,
    interp_dim::Union{String, Number},
    ;
    coord_new = nothing,
    data = nothing,
    data_func::Union{Function, Nothing} = nothing,
    interp_method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    interp_kwargs::Dict = Dict{Symbol, Any}(),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::DCIKT = default_conservative_interp_kwargs,
    weight = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    A::Union{Nothing, AbstractArray} = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    Af::Union{Nothing, AbstractArray} = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    use_svectors::Bool = false,
)
    vardata = isa(var, String) ? data[var] : var
    if ~isnothing(data)
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    else
        interp_dim = isa(interp_dim, String) ? get_dim_num(interp_dim, vardata) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    end

    # evaluate interp_dim_in_is_full_array based on the size of the input... interp_dim_in is full array false is much faster cause dont have to double loop in vectorization...

    if isnothing(weight) # if weight is a number, we assume it's a scalar and we don't need to do anything
        return interp_along_dim(
            vardata,
            interp_dim,
            coord_in;
            interp_dim_out = coord_new,
            data = data,
            data_func = data_func,
            interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
            use_svectors = use_svectors,
        )

    else

        return interp_along_dim(
            vardata .* weight,
            interp_dim,
            coord_in;
            interp_dim_out = coord_new,
            data = data,
            data_func = data_func,
            interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
            use_svectors = use_svectors,
        ) ./ interp_along_dim(
            weight,
            interp_dim,
            coord_in;
            interp_dim_out = coord_new,
            data = data,
            data_func = data_func,
            interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
            use_svectors = use_svectors,
        )

    end
end


"""
    get_data_new_z_t

Take data from our base setup, interpolate it to new z, then create time splines based on the t we have....
to vectorize properly over z_new, it should be the same shape as vardata+vardata_g
"""
function get_data_new_z_t(
    var,
    z_new,
    z_dim::Union{String, Number},
    time_dim::Union{String, Number},
    flight_number::Int;
    thermo_params::TDPS,
    varg = nothing,
    z_old = nothing,
    t_old = nothing,
    data = nothing,
    initial_condition::Bool = false,
    assume_monotonic::Bool = false,
    interp_method::Type{<:AbstractInterpolationMethod} = FastLinear1DInterpolationMethod,
    Spline1D_interp_kwargs::Dict = Dict{Symbol, Any}(),
    pchip_interp_kwargs::Dict = Dict{Symbol, Any}(),
    ground_indices = :end,
    conservative_interp::Bool = false,
    conservative_interp_kwargs::DCIKT = default_conservative_interp_kwargs,
    weight = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    weightg = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    return_before_interp::Bool = false,
    A::Union{Nothing, AbstractArray} = nothing,
    Af::Union{Nothing, AbstractArray} = nothing,
    use_svectors::Bool = false,
)

    if conservative_interp && isnothing(A)
        if interp_method ∈ (DierckxSpline1DInterpolationMethod, FastLinear1DInterpolationMethod)
            A = get_conservative_A(z_new; method = interp_method, Spline1D_interp_kwargs...)
            Af = LinearAlgebra.factorize(A)
        elseif interp_method <: AbstractPCHIPInterpolationMethod
            A = get_conservative_A(z_new; method = interp_method, pchip_interp_kwargs...)
            Af = LinearAlgebra.factorize(A)
        else
            error("unsupported interpolation method $(interp_method) for conservative regridding")
        end
    end


    # do the same processing to the weights and shortcircuit the return w/ return_before_interp
    if !isnothing(weight)
        if !isnothing(weightg)
            weight = get_data_new_z_t(
                weight,
                z_new,
                z_dim,
                time_dim,
                flight_number;
                weight = nothing,
                weightg = nothing,
                varg = weightg,
                return_before_interp = true,
                thermo_params = thermo_params,
                z_old = z_old,
                t_old = t_old,
                data = data,
                initial_condition = initial_condition,
                assume_monotonic = assume_monotonic,
                interp_method = interp_method,
                Spline1D_interp_kwargs = Spline1D_interp_kwargs,
                pchip_interp_kwargs = pchip_interp_kwargs,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                A = A,
                Af = Af,
                use_svectors = use_svectors,
            )
        else
            weight = get_data_new_z_t(
                weight,
                z_new,
                z_dim,
                time_dim,
                flight_number;
                weight = nothing,
                weightg = nothing,
                varg = nothing,
                return_before_interp = true,
                thermo_params = thermo_params,
                z_old = z_old,
                t_old = t_old,
                data = data,
                initial_condition = initial_condition,
                assume_monotonic = assume_monotonic,
                interp_method = interp_method,
                Spline1D_interp_kwargs = Spline1D_interp_kwargs,
                pchip_interp_kwargs = pchip_interp_kwargs,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                A = A,
                Af = Af,
                use_svectors = use_svectors,
            )
        end
    end

    # get the data and dimensions we're working on,
    vardata = isa(var, String) ? data[var] : var
    if ~isnothing(varg)
        vardatag = isa(varg, String) ? data[varg] : varg
    end

    # combine air and ground data
    if ~isnothing(data)
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    else
        z_dim_num = isa(z_dim, String) ? get_dim_num(z_dim, vardata) : z_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
        time_dim_num = isa(time_dim, String) ? get_dim_num(time_dim, vardata) : time_dim
    end


    if isnothing(z_old)
        z_old = z_from_data(data; thermo_params) # uses ground value to create the old z, pads bottom w/ 0
    end
    if isnothing(t_old)
        t_old = data["tsec"][:] # check this unit was right in the files (may need to make sure it's subtracting out the first timestep so starts at 0) -- do we need to align this on a dimension?
    end

    # t_base = Dates.DateTime(string(data["bdate"][:]), Dates.DateFormat("yymmdd")) + Dates.Year(2000) # the base Date (using bdate not nbdate cause nbdate seems to have  bug in flight 9 (extra 0 in month spot))
    # t      = t_base .+ Dates.Second.(t_old) # the actual dates
    # summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
    # SOCRATES_summary = NC.Dataset(summary_file,"r")
    # flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
    # initial_time = SOCRATES_summary["reference_time"][flight_ind] - Dates.Hour(12) # change to select by flight number...
    # initial_ind = argmin(abs.((t.-initial_time))) # find the index of the initial time
    initial_ind = get_initial_ind(data, flight_number, t_old = t_old)

    ### SHOULD WE INTERPOLATE TO THE EXACT TIME RATHER THAN CLOSEST TIMES? IDK... would need to be done before creating vertical splines... 

    if ~isnothing(varg)
        # here we also are gonna need to check where things get inserted in case they are not in order...
        if !assume_monotonic # use data to figure out how and where to do insertions...
        # we need some way to get the local dimension from just a variable
        else
            vardata = combine_air_and_ground_data(vardata, vardatag, z_dim_num; insert_location = ground_indices) # append ground data with 0 as z bottom, loses labeling now though  (this puts a lot of faith im these 2 vars being the same size of having labels which we can't guarantee, no?)
        end
    end

    # select only the initial condition timestep, but keep that dimension around w/ []
    # note -- if not the initial condition, we should still only return initial condition to reference timestep no? (or i guess at least just from the initial condition to the end of the data we have...)
    if initial_condition # the timestep that is closest to the one we are supposed to force with (should be reference time - 12 hours)
        vardata = selectdim(vardata, time_dim_num, [initial_ind])
        z_old = selectdim(z_old, time_dim_num, [initial_ind])
    else # not init condition, so we'll truncate from init condition to end along time dimension
        vardata = selectdim(vardata, time_dim_num, initial_ind:length(t_old))
        z_old = selectdim(z_old, time_dim_num, initial_ind:length(t_old))
    end

    # reverse the z so it goes from ground to top) and matches the new grid we defined..
    z_old = reverse(z_old; dims = z_dim_num) # not reverse!() bc goes SubArray --> Array
    vardata = reverse(vardata; dims = z_dim_num) # not reverse!() bc goes SubArray --> Array

    if return_before_interp
        return vardata
    end



    # interpolate to new z
    if interp_method ∈ (DierckxSpline1DInterpolationMethod, FastLinear1DInterpolationMethod)
        vardata = var_to_new_coord(
            vardata,
            z_old,
            z_dim_num;
            coord_new = z_new,
            data = data,
            interp_method = interp_method,
            interp_kwargs = Spline1D_interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            weight = weight,
            A = A,
            Af = Af,
            use_svectors = false,  # false for internal (z interpolation)
        )
    elseif interp_method <: AbstractPCHIPInterpolationMethod
        vardata = var_to_new_coord(
            vardata,
            z_old,
            z_dim_num;
            coord_new = z_new,
            data = data,
            interp_method = interp_method,
            interp_kwargs = pchip_interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            weight = weight,
            A = A,
            Af = Af,
            use_svectors = false, # false for internal (z interpolation)
        )
    else
        error("unsupported interpolation method")
    end

    if initial_condition # no need to push further here since is init condition (maybe change later to return both?)
        return vardata
    end
    # create new time splines

    # i thnk here should be t_old[initial_ind:end] -- we want to keep only from initial condition timestep
    # when creating the splines, should we use t_old[initial_ind:end] .- t_old[initial_ind]? since the time in the model callling will always start @ t=0

    vardata = var_to_new_coord(
        vardata,
        t_old[initial_ind:end] .- t_old[initial_ind],
        time_dim_num;
        coord_new = nothing,
        data = data,
        interp_method = FastLinear1DInterpolationMethod, # in time, we're gonna stick to linear interpolation for now... this one maybe can be pchip since it's all within the data bounds? idk... i was getting w=0 using pchip... possibly because
        interp_kwargs = Spline1D_interp_kwargs,
        conservative_interp = false, # no need for conservation in time? maybe?
        use_svectors = use_svectors, # potentially true for time interpolation
    )

    return vardata
end



function drop_lat_lon(vardata; data = nothing, dims = nothing)
    # gotta do this at the end after we've made full use of our dimnames since our data is unlabeled
    # in general we could cheat in the future since the socrates order seems to always be lon lat lev time regardless of which dims exist...
    # i guess i also don't know what happens if we've turned the time dim into just a fcn -- in principle it's last so that shouldnt hurt

    if isnothing(dims)
        if isnothing(data)
            throw(ArgumentError("You must provide `data` if `dims` is not specified"))
        end
        dims = tuple(collect(findfirst(x -> x == dim, NC.dimnames(data["T"])) for dim in ["lat", "lon"])...) # base off temperature for now
    end

    vardata = dropdims(vardata, dims = dims)
    return vardata
end

"""
    insert_dims

add in new dimensions at position (never seemed to use , could get rid of this?)
"""
function insert_dims(data, ind; new_dim_sizes = [-1])
    sz_data = collect(size(data)) # array
    splice!(sz_data, ind, [new_dim_sizes..., sz_data[ind]]) # do this way cause splice itself deletes the current value so preserve it
    data = reshape(data, sz_data...) # reshape
end

"""
Currently this is setup to just assume saturation w/ Tg 
However, this doesn't match Atlas's simulations so maybe we'll adjust this to just be an adiabatic adjustment to the lowest datapoint we do have?
It seems that Atlas's simulations have a slight kink at the lowest level, but otherwise are ≈ constant down to the surface if that's outside the forcing range and just interpolated if it's inside the forcing range...
    - in that case, we should be able to just use interpolate_1d (used to be pyinterp) because Spline1D default bc is nearest outside the range
"""
# function calc_qg(Tg,pg; thermo_params)
function calc_qg_extrapolate_pq(pg, p, q; interp_method = FastLinear1DInterpolation)
    # pvg           = TD.saturation_vapor_pressure.(thermo_params, Tg, TD.Liquid())
    # molmass_ratio = TDP.molmass_ratio(thermo_params)
    # qg            = (1 / molmass_ratio) .* pvg ./ (pg .- pvg) #Total water mixing ratio at surface , assuming saturation [ add source ]

    # not sure if this should be linear in p or logarithmic (linear in z), gonna do linear in p
    qg = interpolate_1d(pg, p, q, interp_method; bc = ExtrapolateBoundaryCondition()) # default to spline1d for extrapolation as pchip may not be the most reliable outside bounds of data.
    # qg = interpolate_1d(log.(pg), log.(p), q, interp_method; bc="extrapolate") # default to spline1d for extrapolation as pchip may not be the most reliable outside bounds of data.

    return qg
end

function calc_qg_from_pgTg(pg, Tg, thermo_params::TDPS)
    ρg = TD.air_density(thermo_params, Tg, pg) # original T?
    qg = TD.q_vap_saturation_generic(thermo_params, Tg, ρg, TD.Liquid()) # surface specific humidity over liquid
    return qg
end

"""
Retrieve the starting index for our LES simulations (12 hours before reference file) in the forcing data
"""
function get_initial_ind(data, flight_number::Int; t_old = nothing)
    if isnothing(t_old)
        t_old = data["tsec"][:] # check this unit was right in the files (may need to make sure it's subtracting out the first timestep so starts at 0) -- do we need to align this on a dimension?
    end

    t_base = Dates.DateTime(string(data["bdate"][:]), Dates.DateFormat("yymmdd")) + Dates.Year(2000) # the base Date (using bdate not nbdate cause nbdate seems to have  bug in flight 9 (extra 0 in month spot))
    t = t_base .+ Dates.Second.(t_old) # the actual dates

    summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
    SOCRATES_summary = NC.Dataset(summary_file, "r")

    flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
    initial_time = SOCRATES_summary["reference_time"][flight_ind] - Dates.Hour(12) # change to select by flight number...
    initial_ind = argmin(abs.((t .- initial_time))) # find the index of the initial time
    return initial_ind
end

function squeeze(x)
    return dropdims(x, dims = (findall(size(x) .== 1)...,))
end

set_property(x::NamedTuple, field::Symbol, value) = merge(x, (field => value,))  # edit the value in a named tuple
