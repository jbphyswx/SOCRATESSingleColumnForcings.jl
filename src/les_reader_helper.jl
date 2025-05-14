#=

Read LES output file instead of the forcing input to set things

=#




"""
This function is for reading things from LES since they don't have the same time setup and already have a z (which should match our desired output)
LES data only has dimensions (time, z), no lev etc, but we keepthe syntax the same to make switching bacc n forth easier
"""
function get_data_new_z_t_LES(
    var,
    z_new, # you'd htink this would be the same as our desired output z, but RF09 was run on a taller grid for some reason
    z_dim::Union{Int, String},
    time_dim::Union{Int, String},
    flight_number::Int,
    forcing_type::Symbol;
    varg = nothing,
    z_old = nothing,
    t_old = nothing,
    data = nothing,
    initial_condition::Bool = false,
    assume_monotonic::Bool = false,
    interp_method::Symbol = :Spline1D,
    Spline1D_interp_kwargs::Dict = Dict{Symbol, Any}(:bc => "extrapolate"), # default to extrapolate bc RF09 was run on a different grid for some reason, we're not necessarily guaranteed to have the same z
    pchip_interp_kwargs::Dict = Dict{Symbol, Any}(:bc => "extrapolate"),
    ground_indices = :end,
    conservative_interp::Bool = false,
    conservative_interp_kwargs::DCIKT = default_conservative_interp_kwargs,
    weight = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    weightg = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    return_before_interp::Bool = false,
    A::Union{Array, Nothing} = nothing,
    Af::Union{Array, Nothing} = nothing,
)

    if conservative_interp && isnothing(A)
        LinearAlgebra.BLAS.set_num_threads(1) # if you're on HPC this is essential to A\b not slowing down by 5 orders of magnitude from 1ms to 100s
        if interp_method ∈ [:Spline1D, :Dierckx]
            A = get_conservative_A(z_new; method = interp_method, Spline1D_interp_kwargs...)
        elseif interp_method ∈ [:pchip_smooth_derivative, :pchip_smooth]
            A = get_conservative_A(z_new; method = interp_method, pchip_interp_kwargs...)
        else
            error("unsupported interpolation method")
        end
        Af = LinearAlgebra.factorize(A) # factorize the matrix for faster solving
    end

    # Preprocess the weights by running the function with the weights as the data but shortcircuiting the response usnig return_before_interp
    if !isnothing(weight)
        if !isnothing(weightg)
            weight = get_data_new_z_t_LES(
                weight,
                z_new,
                z_dim,
                time_dim,
                flight_number,
                forcing_type;
                weight = nothing,
                weightg = nothing,
                varg = weightg,
                return_before_interp = true,
                z_old = z_old,
                t_old = t_old,
                data = data,
                initial_condition = initial_condition,
                assume_monotonic = assume_monotonic,
                interp_method = interp_method,
                Spline1D_interp_kwargs = Spline1D_interp_kwargs,
                pchip_interp_kwargs = pchip_interp_kwargs,
                ground_indices = ground_indices,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                A = A,
                Af = Af,
            )
        else
            weight = get_data_new_z_t_LES(
                weight,
                z_new,
                z_dim,
                time_dim,
                flight_number,
                forcing_type;
                weight = nothing,
                weightg = nothing,
                varg = nothing,
                return_before_interp = true,
                z_old = z_old,
                t_old = t_old,
                data = data,
                initial_condition = initial_condition,
                assume_monotonic = assume_monotonic,
                interp_method = interp_method,
                Spline1D_interp_kwargs = Spline1D_interp_kwargs,
                pchip_interp_kwargs = pchip_interp_kwargs,
                ground_indices = ground_indices,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                A = A,
                Af = Af,
            )
        end
    end


    if isa(var, String) && isnothing(data)
        data = open_atlas_les_output(flight_number, forcing_type)[forcing_type]
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
        z_old = data["z"][:] # z is already in LES data
    end
    if isnothing(t_old) # is in unis of days in LES files
        t_old = data["time"][:] # check this unit was right in the files (may need to make sure it's subtracting out the first timestep so starts at 0) -- do we need to align this on a dimension?
        t_old = (t_old .- t_old[1]) * 24 * 3600 # convert to seconds
    end

    initial_ind = 1 # initial ind is just 0 since we're using th LES output
    if isa(vardata, NC.CFVariable) # breaks things downstream if it is and we aren't calling combine air and ground data like we did for the others..
        vardata = vardata[:]
    end

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
    else # not init condition, so we'll truncate from init condition to end along time dimension
        vardata = selectdim(vardata, time_dim_num, initial_ind:length(t_old))
    end

    # reverse the z so it goes from ground to top) and matches the new grid we defined.. (is this ncessary)
    # z_old = reverse(z_old; dims = z_dim_num)
    # vardata = reverse(vardata; dims = z_dim_num)

    if return_before_interp # short circuit and just return the preprocessed vardata
        return vardata
    end

    # interpolate to new z (though it should alreayd be on the new z, though I guess you could do some smoothing/rounding etc)

    if interp_method ∈ [:Spline1D, :Dierckx]
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
        ) # extrapolate bc it's not a guarantee that our new z  will contain our desired z (mostly bc RF09 was run on a different grid for some reason) 
    elseif interp_method ∈ [:pchip_smooth_derivative, :pchip_smooth]
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
        t_old[initial_ind:end] .- t_old[initial_ind], # techincally this should be uncessesary for LES
        time_dim_num;
        coord_new = nothing,
        data = data,
        interp_method = :Spline1D, # in time, we're gonna stick to linear interpolation for now... this one maybe can be pchip since it's all within the data bounds? idk... i was getting w=0 using pchip... possibly because
        interp_kwargs = Spline1D_interp_kwargs,
        conservative_interp = false, # not need for conservative interp in time... maybe? it's philosophical and depends on the variable... does the time-averaged or instantaneous value matter more... However, we are just creating time splines and the timestep is up to your model... your model thus has to handle any time averaging etc...
    )

    return vardata

end
