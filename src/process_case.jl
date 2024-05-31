
"""
    process_case(
        flight_number::Int;
        obs_or_ERA5 = "Obs"::Union{String,Symbol},
        new_z::Union{Nothing,AbstractArray}=nothing,
        initial_condition::Bool=false,
        thermo_params,
        surface=nothing
    )

Processes the flight data for that case
If new_z is not specified, will use Rachel Atlas default grid
If initial_condition is false, returns the data at time 0, if true, returns full spline interpolations over the time dimension

Some things are always forced by ERA5, otherwise we pull from whatever obs_or_ERA5 specifies... (list which is which?)
"""
function process_case(
    flight_number::Int;
    obs_or_ERA5 = "Obs"::Union{String, Symbol},
    new_z::Union{Nothing, AbstractArray, NamedTuple} = nothing,
    initial_condition::Bool = false,
    thermo_params,
    surface = nothing,
    use_LES_output_for_z = false,
    return_old_z = false,
)

    # initial conditions
    data = open_atlas_les_input(flight_number)

    # get the new grid we want from the atlas les grid file
    return_values = (:dTdt_hadv, :H_nudge, :dqtdt_hadv, :qt_nudge, :subsidence, :u_nudge, :v_nudge, :ug_nudge, :vg_nudge, :dTdt_rad)
    if isnothing(new_z)
        new_z = data[:grid_data]
        # named tuple repeated for each in return values
        new_z = NamedTuple{return_values}((new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z))
    elseif isa(new_z, AbstractArray)
        new_z = NamedTuple{return_values}((new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z, new_z))
    end

    # resolve forcing keyword
    if obs_or_ERA5 ∈ ["Obs", :obs_data]
        forcing = :obs_data
        # data = data[(:obs_data, :ERA5_data)]
        data = data[(:obs_data,)] # drop ERA if we're doing obs cause we don't need it (11 has no obs either i think)        
    elseif obs_or_ERA5 ∈ ["ERA5", :ERA5_data]
        forcing = :ERA5_data
        data = data[(:ERA5_data,)] # drop obs if we're doing era5 cause we don't need it (11 has no obs either i think)        
    else
        if obs_or_ERA5 ∉ [:obs_data, :ERA5_data, "Obs", "ERA5"]
            error("obs_or_ERA5 must be either \"Obs\"/[:obs_data] or \"ERA5\"/[:ERA5_data]")
        else
            forcing = obs_or_ERA5
        end
    end

    # specify our action dimensions
    z_dim_num = get_dim_num("lev", data[forcing]["T"])
    time_dim_num = get_dim_num("time", data[forcing]["T"])

    # We use map here a lot to map over all forcings in our data, but that could maybe be deprecated if we stop carrying it around as a named tuple

    ## For surface values, adjust to the time period of interest and return only the surface values, these are called in TC.jl
    if ~isnothing(surface)
        initial_ind = get_initial_ind(data[forcing], flight_number) # should the reference be the first timestep? or what is the reference meant to be...
        summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
        SOCRATES_summary = NC.Dataset(summary_file, "r")
        flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
        Tg_offset = SOCRATES_summary[:deltaT][flight_ind]  # I think this is backwards in table 2 in the paper... is really T_2m - SST as in Section 3
        if surface ∈ ["reference_state", "reference", "ref"]  # we just want the surface reference state and we'll just return that
            Tg = data[forcing]["Tg"][:][initial_ind] + Tg_offset # might have to drop lon,lat dims or sum
            pg = data[forcing]["Ps"][:][initial_ind]

            if Tg_offset < 0  # SST/T_orig > Tg, assume SST sets qg at ground level and serves as a source
                qg = calc_qg_from_pgTg(pg, Tg, thermo_params)
            else # SST/T_orig < Tg, stable boundary layer, assume latent air RH controls things... not SST
                p = vec(data[forcing]["lev"])[:]
                q = vec(selectdim(data[forcing]["q"], time_dim_num, initial_ind))[:] # select our q value subset along the time dimension
                q = q ./ (1 .+ q) # mixing ratio to specific humidity
                qg = calc_qg_extrapolate_pq([pg], p, q)
                qg = collect(qg)[]  # Thermodynamics 0.10.2 returns a tuple rather than scalar, so this can collapse to scalar in either 0.10.1<= or 0.10.2>=
            end
            return TD.PhaseEquil_pTq(thermo_params, pg, Tg, qg)
        elseif surface ∈ ["surface_conditions", "conditions", "cond"]
            pg = vec(data[forcing]["Ps"])[:][initial_ind:end]
            Tg_orig = vec(data[forcing]["Tg"])[:][initial_ind:end] # SST
            Tg = Tg_orig .+ Tg_offset # might have to drop lon,lat dims or sum
            # before we were extrapolating to surface to get qg which looks ok at first but after ts 1, the LES diverges, maybe it's supposed to be going towards Tg SST
            if Tg_offset < 0  # SST > Tg, assume SST sets qg at ground level going forward and serves as a source (so use full Tg not Tg_orig)
                qg = calc_qg_from_pgTg.(pg, Tg, thermo_params)  # surface specific humidity over liquid
            else # SST < Tg, so stable boundary layer, so still set qg to the exigent air RH value not the SST saturation so that it doesnt serve as a source of moisture (even if SST saturation value is greater) (Not sure why boundary layer parameterization doesn't just handle it but we were getting huge moisture spikes at sfc from SST when RH was low, maybe it's a subsidence thing)
                p = vec(data[forcing]["lev"])[:]
                q = selectdim(data[forcing]["q"], time_dim_num, initial_ind:size(data[forcing]["q"], time_dim_num)) # select our q value subset along the time dimension
                q = vec.(collect(eachslice(q, dims = time_dim_num))) # turn our q from [lon, lat, lev, time] to a list of vectors along [lon,lat,lev] to match p
                q = map(mr -> mr ./ (1.0 .+ mr), q) # mixing ratio to specific humidity for each vector we created in q
                qg = map((pg, q) -> calc_qg_extrapolate_pq([pg], p, q)[1], pg, q) # map the function to get out qg for each time step
            end

            tg = data[forcing]["tsec"][initial_ind:end] # get the time array
            tg = tg .- tg[1] # i think we need this to get the initial time to be 0, so the interpolation works
            # in this interpolation, tsec has to be adjusted to our offsets no? or we clip e.g. pg to be pg[initial_ind:end], tsec would also need to be adjusted no?, subtract the value at initial_ind i guess...
            return (;
                pg = t -> pyinterp([t], tg, pg; method=:Spline1D)[1], # in time always use spline1d rn...
                Tg = t -> pyinterp([t], tg, Tg; method=:Spline1D)[1],
                qg = t -> pyinterp([t], tg, qg; method=:Spline1D)[1],
            ) # would use ref and broadcast but doesnt convert back to array
        else
            error(
                "if surface is not set to nothing (i.e you want a surface value, it must be either \"reference_state\" or \"surface_conditions\"",
            )
        end
    end

    # p,T,q | need these from both ERA and obs to construct ts, tsg
    p = map(x -> x["lev"], data)
    p = map(x -> align_along_dimension(x, z_dim_num), p) # align on dimension 3 lev (hopefully order was already correct (testing, can't recall if theres anywhere i didn't want this done...)
    pg = map(x -> x["Ps"], data)
    T = map(x -> x["T"], data)
    Tg = map(x -> x["Tg"], data)
    q = map(x -> x["q"], data)
    q = map(mr -> mr ./ (1 .+ mr), q) # mixing ratio to specific humidity for each forcing

    #For not surface though, we need to add in the ΔT from the summary table in the Atlas paper (to tsg above?)
    summary_file = joinpath(dirname(@__DIR__), "Data", "SOCRATES_summary.nc")
    SOCRATES_summary = NC.Dataset(summary_file, "r")
    flight_ind = findfirst(SOCRATES_summary["flight_number"][:] .== flight_number)
    Tg_orig = map(x -> x, Tg) # copy (can we just use copy()?)
    Tg_offset = SOCRATES_summary[:deltaT][flight_ind]  # I think this is backwards in table 2 in the paper... is really T_2m - SST as in Section 3
    Tg = map(Tg -> Tg .+ Tg_offset, Tg)

    if Tg_offset < 0  # SST > Tg, assume Tg limits moisture below last known point and just extrapolate
        base_calc_qg = (pg, p, q) -> calc_qg_extrapolate_pq([pg], p, q[:])[1] # pg->[pg] for pyinterp and [1] for just the value out
        qg = map(
            (pg, p, q) ->
                base_calc_qg.(
                    pg,
                    Ref(p[:]),
                    align_along_dimension(vec.(collect(eachslice(q; dims = time_dim_num))), z_dim_num),
                ),
            pg,
            p,
            q,
        ) # iterate over forcings, the pg value, the p value is a fixed array, for the q value we take our slices in z and align them along the time dimension to match the shape of pg for calc_qg broadcasting, 

    else # Tg > SST/T_orig, assume SST sets moisture  below last known point (t > t_orig), t_orig = Tg_q_sfc
        # Tg_q_sfc = map((x,y) -> min.(x, y), Tg, Tg_orig) # min of Tg and T (Moisture comes from evaporation, if Tg < SST, SST limits moisture, if Tg > SST, Tg limits moisture) This is important bc it's how we extrapolate from the surface...
        base_calc_qg = (pg, Tg) -> calc_qg_from_pgTg(pg, Tg, thermo_params) # pg->[pg] for pyinterp and [1] for just the value out 
        qg = map(
            (pg, Tg) -> base_calc_qg.(pg, Tg),
            pg,
            Tg_orig, ### CHECK WHETHER WE WANT TO BE USING BASE TG OR OFFSET TG ### (OR MAYBE DIFF AT T=0 VS AFTERWARDS? WE DON'T ACTUALLY EVER SET TG (Tg_orig was wayyy too high on RF09, Tg is too high on RF10)
        ) # TESTINGGGG (if p_input < pg, there may be a sfc discontinuity, otherwise we should get cleaner extrapolation...
    end

    # Set up thermodynamic states for easier use (for both forcing and ERA -- note ERA subsidence for example depends on density which relies on T,p,q so need both even if forcing is :obs_data)
    tsg = map((pg, Tg, qg) -> TD.PhaseEquil_pTq.(thermo_params, pg, Tg, qg), pg, Tg, qg)

    # In principle at initiation they assume domination by liquid, so T_l,i ≈ T_l. We are given in the forcing files T_L = T - L q_c/ c_p. Then, θ_L,I ≈ θ_L = θ - θ\T L q_c/ c_p = (θ/T) ( T - L q_c/ c_p) = (θ/T) T_L thus θ_L = T_L (p_0/p)^k, the same as just calculating the dry potential temperature subsitututing T_L
    θ = map((T, p,) -> TD.dry_pottemp_given_pressure.(thermo_params, T, p), T, p) # should be same as # θ = map((T, p, q) -> TD.dry_pottemp_given_pressure.(thermo_params, T, p, TD.PhasePartition.(q)), T, p, q)

    ts = map((p, θ, q) -> TD.PhaseEquil_pθq.(thermo_params, p, θ, q), p, θ, q)

    # get indices where the ground would get inserted in... (and convert to same shape/dims.)
    ground_indices = map((ts) -> get_ground_insertion_indices(ts, tsg[forcing], z_dim_num; thermo_params, data = data), ts) # 

    q_full = combine_air_and_ground_data(q[forcing], qg[forcing], z_dim_num; insert_location = ground_indices[forcing]) # full q_array/qt_nudge -- is used from the chosen forcing dataset
    qt_nudge = q_full


    ts_full = map(
        (ts, ground_indices) ->
            combine_air_and_ground_data(ts, tsg[forcing], z_dim_num; insert_location = ground_indices),
        ts,
        ground_indices,
    ) # is it era5? their output les plots sure don't look it...

    # SHOULD THESE NOT USE THE GROUND_INDICES INSERT LOCATIONS???? and/or ts_full

    # old_z  => Precompute old z coordinate (precompute to save us some trouble later (get_data_new_z_t func can self-calculate it but it's redundant to keep calculating z)
    if !use_LES_output_for_z
        z_old = map((ts, data) -> lev_to_z(ts, tsg[forcing]; thermo_params, data = data), ts, data) # should this be tsg[:ERA5_data] cause surface is always ERA5 (is it?) 
    else
        z_old = map(
            (ts, data, ground_indices) -> lev_to_z_from_LES_output(
                ts,
                tsg[forcing];
                thermo_params,
                data = data,
                flight_number = flight_number,
                forcing_type = forcing,
                ground_indices = ground_indices,
            ),
            ts,
            data,
            ground_indices,
        ) # should this be tsg[:ERA5_data] cause surface is always ERA5 (is it?) [ + Testing getting z from the forcing data ]
    end

    if return_old_z
        return z_old
    end


    # ω (subsidence) # always forced by era5
    ρ = TD.air_density.(thermo_params, ts[forcing])
    ρg = TD.air_density.(thermo_params, tsg[forcing])
    ρ = combine_air_and_ground_data(ρ, ρg, z_dim_num; insert_location = ground_indices[forcing])
    ω = data[forcing]["omega"]
    dpdt_g = data[forcing]["Ptend"]
    dpdt_g = add_dim(dpdt_g, z_dim_num) # should be lon lat lev time (hopefully order was already correct)
    ω = combine_air_and_ground_data(ω, dpdt_g, z_dim_num; insert_location = ground_indices[forcing])

    p_grid = map((p,) -> add_dim(align_along_dimension(p, z_dim_num), time_dim_num), p) # align on dimension 3 lev, add time dimensino 4
    p_full = map(
        (p_grid, pg, ground_indices) ->
            combine_air_and_ground_data(p_grid, pg[:], z_dim_num; insert_location = ground_indices),
        p_grid,
        pg,
        ground_indices,
    ) # align p along lev, add lev_dim to pg, stack (hope braodcasting works to fill it out...)

    # this i think was wrong -- maybe get atlas to confirm what version of the sigmoid she used.
    # c    = 100
    # a    = -2 + exp(250/c)
    # f_p  = @. 2(a+1) / (a+exp(p_full/c)) - 1

    # squeeze = (x) -> dropdims(x, dims = (findall(size(x) .== 1)...,))

    # this was derived for a scaled sigmoid passing through (ps,1), (25000 Pa, 0)
    L = 2.2 # maximum value (shape parameter)
    a = -L / 2
    p0 = 250.0 * 100
    p1 = add_dim(pg[forcing][:], z_dim_num)
    (p1, y1) = (p1, 1) # use our ground pressure (can't use pfull cause would need to pull)
    (p2, y2) = (p0, 0)
    k = @. log((L / 2 + 1) / (L / 2 - 1)) / (p1 - p0) # ps is an array so we have an array of ks
    f_p = @. cos(π / 2 * (p1 - p_full[forcing]) / (p1 - p0)) * (p_full[forcing] >= p0) # atlas email
    f_p_alt = @. (a + L / (1 + exp(-k * (p_full[forcing] - p0)))) * (p_full[forcing] >= p0) # my original



    grav = TDP.grav(thermo_params)
    subsidence = -(ω .- (dpdt_g .* f_p)) ./ (ρ .* grav)


    # maybe it's an upwinding thing? what if we turn off ascent within 3 indices of either edge
    # subsidence_buffer = 0
    # nz = size(subsidence, z_dim_num)
    # valid_subsidence = align_along_dimension(1:nz, z_dim_num)
    # valid_subsidence = (valid_subsidence .> subsidence_buffer) .& (valid_subsidence .<= (nz - subsidence_buffer))
    # subsidence = subsidence .* (subsidence .>0) + subsidence .* (subsidence .<0) .* valid_subsidence # remove any ascent near boundaries to test if that is more stable...

    # subsidence = max.(subsidence,0)  # stability test (unstable) (all ascent)
    # subsidence = min.(subsidence, 0)  # stability test (stable!) ( all subsidence ) # testing getting rid of this to see if the stability problem had been fixed elsewhere (it hasn't lol, tested)
    # m = .007 # almost stable .005 | unstable .007
    # subsidence[subsidence .> m] .= m # stability test -- can we tolerate any ascent?


    # u, v # always forced by ERA5
    u = data[forcing]["u"]
    v = data[forcing]["v"]
    u = combine_air_and_ground_data(u, FT(0), z_dim_num; insert_location = ground_indices[forcing])
    v = combine_air_and_ground_data(v, FT(0), z_dim_num; insert_location = ground_indices[forcing])
    u_nudge, v_nudge = u, v

    #=
    not sure what
        "forced by ERA5-derived geostrophic winds and nudged toward the ERA5 horizontal winds with a nudging
        timescale of 1 hr for the ERA5-based simulation and 20 min for the Obs-based simulation."
    means so we're just using u,v not ug,vg for now...
    What does it mean to have both the forcing and the nudging? are they not overlapping/competing?
    =#

    # u_g, v_g # always forced by ERA5
    ug = data[forcing]["ug"]
    vg = data[forcing]["vg"]
    ug = combine_air_and_ground_data(ug, FT(0), z_dim_num; insert_location = ground_indices[forcing])
    vg = combine_air_and_ground_data(vg, FT(0), z_dim_num; insert_location = ground_indices[forcing])
    ug_nudge, vg_nudge = ug, vg

    # H (nudge)
    θ_liq_ice = TD.liquid_ice_pottemp.(thermo_params, ts_full[forcing])
    H_nudge = θ_liq_ice
   

    # dTdt_hadv (i think these are supposed to be always ERA5)
    dTdt_hadv = data[forcing]["divT"]
    dTdt_hadv = combine_air_and_ground_data(dTdt_hadv, FT(0), z_dim_num; insert_location = ground_indices[forcing])

    # dqdt_hadv (i think these are supposed to be always ERA5)
    dqtdt_hadv = data[forcing]["divq"]
    dqtdt_hadv = combine_air_and_ground_data(dqtdt_hadv, FT(0), z_dim_num; insert_location = ground_indices[forcing])
    # expand everything to new z grid and make t operation (initial condition or time splines)
    # dTdt_hadv  = get_data_new_z_t(dTdt_hadv , new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)

    # ======================================================================================================================== #
    LES_data = open_atlas_les_output(flight_number)[forcing]
    dTdt_rad = LES_data["RADQR"] ./ (24 * 3600) # convert to K/s from K/day (dividing drops nc dim data tho but converts it to array from ncvariable)

    z_dim_num_LES = get_dim_num("z", LES_data["RADQR"]) # assume it's same for all LES 2D vars?
    time_dim_num_LES = get_dim_num("time", LES_data["RADQR"])
    
    dTdt_rad = get_data_new_z_t_LES(
        dTdt_rad,
        new_z[:dTdt_rad],
        z_dim_num_LES,
        time_dim_num_LES,
        flight_number;
        z_old = nothing,
        data = LES_data,
        thermo_params,
        initial_condition,
    )[:]

    # ======================================================================================================================== #

    dTdt_hadv = get_data_new_z_t(
        dTdt_hadv ,#.* (forcing == :ERA5_data),
        new_z[:dTdt_hadv],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:]  # testing not using forcing for obs to see if fixes anything

    H_nudge = get_data_new_z_t(
        H_nudge,
        new_z[:H_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
        # interp_method = :pchip_smooth, # testing if this is the one that breaks things
        interp_method = :Spline1D,
        pchip_interp_kwargs = Dict(
            :f_enhancement_factor => 5, # higher to keep sharp inversions
            :f_p_enhancement_factor => 8),  # not too high to avoid cusps (changed to high to keep model fidelity)
    )[:]

    # dqtdt_hadv = get_data_new_z_t(dqtdt_hadv, new_z, z_dim_num,time_dim_num, flight_number; z_old = z_old[forcing], data=data[forcing], thermo_params,  initial_condition)
    dqtdt_hadv = get_data_new_z_t(
        dqtdt_hadv ,#.* (forcing == :ERA5_data),
        new_z[:dqtdt_hadv],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:] # testing not using forcing for obs to see if fixes anything

    qt_nudge = get_data_new_z_t(
        qt_nudge,
        new_z[:qt_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
        # interp_method = :pchip_smooth,
        interp_method = :Spline1D,
        pchip_interp_kwargs = Dict(
            :f_enhancement_factor => 6, # higher to keep sharp inversions
            :f_p_enhancement_factor => 8),  # not too high to avoid cusps
    )[:]

    subsidence = get_data_new_z_t(
        subsidence,
        new_z[:subsidence],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
        interp_method = :Spline1D, # verifying if this is still unstable
        # interp_method = :pchip_smooth,
        pchip_interp_kwargs = Dict(
            :f_enhancement_factor => 1, # lower for gentle changes and no sharp convergence/divergence, loss of accuracy ok
            :f_p_enhancement_factor => 1),  # lower for gentle changes and no sharp convergence/divergence, loss of accuraacy ok
    )[:]

    u_nudge = get_data_new_z_t(
        u_nudge,
        new_z[:u_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:]

    v_nudge = get_data_new_z_t(
        v_nudge,
        new_z[:v_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:]

    ug_nudge = get_data_new_z_t(
        ug_nudge,
        new_z[:ug_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:]
    vg_nudge = get_data_new_z_t(
        vg_nudge,
        new_z[:vg_nudge],
        z_dim_num,
        time_dim_num,
        flight_number;
        z_old = z_old[forcing],
        data = data[forcing],
        thermo_params,
        initial_condition,
    )[:]

    return (; dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge, ug_nudge, vg_nudge, dTdt_rad)

end


function get_default_new_z(flight_number::Int;)
    data = open_atlas_les_input(flight_number)
    new_z = data[:grid_data]
    return new_z
end

"""
    surface_ref_state(
        flight_number::Int;
        obs_or_ERA5 = "Obs"::Union{String,Symbol},
        thermo_params
    )

In Cases.jl from TurbulenceConvection.jl you're asked to provide a surface reference state so I'm gonna provide one based off just whatever our base profile is
"""
function surface_ref_state(flight_number::Int; obs_or_ERA5 = "Obs"::Union{String, Symbol}, thermo_params)
    return nothing
end
