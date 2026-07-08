# Field-array altitude / ground-insertion helpers (no thermodynamic state carriers).

"""
    ground_insertion_indices(p, pg, concat_dim; data)

Pressure-based ground insertion indices
"""
function ground_insertion_indices(
    p::AbstractArray,
    pg::AbstractArray,
    concat_dim::Union{Int, String};
    data,
)
    # Materialize lazy disk-backed arrays so reshape/cat/mapslices use stable Array methods.
    p_arr = Array(p)
    pg_arr = Array(pg)
    concat_dim = isa(concat_dim, String) ? dim_num(concat_dim, p_arr) : concat_dim

    sz_p = collect(size(p_arr)) # this allocates
    sz_pg = collect(size(pg_arr)) # this sallocates
    if ndims(pg_arr) == ndims(p_arr) - 1
        insert!(sz_pg, concat_dim, 1)
        pg_shaped = reshape(pg_arr, sz_pg...)
    elseif ndims(pg_arr) == ndims(p_arr)
        pg_shaped = pg_arr
    else
        error(
            "ground_insertion_indices: pg ndims=$(ndims(pg_arr)) incompatible with p ndims=$(ndims(p_arr))",
        )
    end

    num_repeat = ones(Int, length(sz_p)) # this allocates
    num_repeatg = ones(Int, length(sz_p)) # this allocates
    for i in eachindex(sz_p)
        i == concat_dim && continue
        sz_a = sz_p[i]
        sz_g = size(pg_shaped, i)
        if sz_a == sz_g
            continue
        elseif sz_a == 1
            num_repeat[i] = sz_g
        elseif sz_g == 1
            num_repeatg[i] = sz_a
        else
            error(
                "ground_insertion_indices: size mismatch in dimension $i (p=$sz_a, pg=$sz_g)",
            )
        end
    end
    p_arr = repeat(p_arr, num_repeat...)
    pg_shaped = repeat(pg_shaped, num_repeatg...)

    combined = cat(p_arr, pg_shaped; dims = concat_dim)

    function mapslice_func(vect)
        vardata = collect(vect[1:(end - 1)])
        filter!(!ismissing, vardata)
        vardatag = vect[end]
        if ismissing(vardatag)
            return length(vardata) + 1
        end
        return searchsortedfirst(vardata, vardatag)
    end

    return mapslices(mapslice_func, combined; dims = [concat_dim])
end

"""
    lev_to_z_column_pTq(p_vec, T_vec, q_vec; thermo_params)

Hypsometric integration on one column; last element of each vector is the ground value.
"""
function lev_to_z_column_pTq(
    p_vec::AbstractVector,
    T_vec::AbstractVector,
    q_vec::AbstractVector;
    thermodynamics_backend =DefaultThermodynamicsBackend(),
)
    FT = promote_type(eltype(p_vec), eltype(T_vec), eltype(q_vec), Float64)
    _R_d = R_d(thermodynamics_backend, FT)
    _grav = grav(thermodynamics_backend, FT)

    L = length(p_vec) - 1
    pg = FT(p_vec[end])
    Tg = FT(T_vec[end])
    qg = FT(q_vec[end])

    n_air = 0
    @inbounds for i in 1:L
        !ismissing(p_vec[i]) && (n_air += 1)
    end

    p_work = Vector{FT}(undef, n_air + 1) # this allocates, bad
    T_work = Vector{FT}(undef, n_air + 1) # this allocates, bad
    q_work = Vector{FT}(undef, n_air + 1) # this allocates, bad

    j = 0
    @inbounds for i in 1:L
        ismissing(p_vec[i]) && continue
        j += 1
        p_work[j] = FT(p_vec[i])
        T_work[j] = FT(T_vec[i])
        q_work[j] = FT(q_vec[i])
    end

    index = searchsortedfirst(view(p_work, 1:j), pg)
    @inbounds for k in j:-1:index
        p_work[k + 1] = p_work[k]
        T_work[k + 1] = T_work[k]
        q_work[k + 1] = q_work[k]
    end
    p_work[index] = pg
    T_work[index] = Tg
    q_work[index] = qg
    n_work = j + 1

    Tvz = Vector{FT}(undef, n_work) # This allocates, bad
    @inbounds for i in 1:n_work
        q_liq, q_ice = equilibrium_condensate(thermodynamics_backend, T_work[i], p_work[i], q_work[i])
        Tvz[i] = virtual_temperature(thermodynamics_backend, T_work[i], q_work[i], q_liq, q_ice)
    end

    dz = Vector{FT}(undef, n_work - 1)
    @inbounds for i in 1:(n_work - 1)
        Tv_bar = (Tvz[i] + Tvz[i + 1]) / 2 # layer-mean virtual temperature (hypsometric integrand)
        dz[i] = (_R_d * Tv_bar / _grav) * log(p_work[i + 1] / p_work[i])
    end

    z = Vector{FT}(undef, n_work) # This allocates, bad
    z[n_work] = zero(FT)
    @inbounds for i in (n_work - 1):-1:1
        z[i] = z[i + 1] + dz[i]
    end
    return z .- z[index] # This allocates, bad
end

"""
    lev_to_z(p, T, q, pg, Tg, qg; data, assume_monotonic=false, ground_indices=nothing)

Field-array altitude from pressure levels and surface fields.
"""
function lev_to_z(
    p::AbstractArray,
    T::AbstractArray,
    q::AbstractArray,
    pg::AbstractArray,
    Tg::AbstractArray,
    qg::AbstractArray;
    data,
    thermodynamics_backend =DefaultThermodynamicsBackend(),
    assume_monotonic::Bool = false,
    ground_indices = nothing,
)
    dimnames = NC.dimnames(data["T"])
    lev_dim_num = findfirst(x -> x == "lev", dimnames)
    ldn = lev_dim_num

    if assume_monotonic
        error("assume_monotonic lev_to_z not implemented for field arrays")
    end

    p_full = combine_air_and_ground_data(p, pg, ldn; data, reshape_ground = true, insert_location = :end)
    T_full = combine_air_and_ground_data(T, Tg, ldn; data, reshape_ground = true, insert_location = :end)
    q_full = combine_air_and_ground_data(q, qg, ldn; data, reshape_ground = true, insert_location = :end)

    pqT = map((pp, TT, qq) -> (Float64(pp), Float64(TT), Float64(qq)), p_full, T_full, q_full)

    return mapslices(col -> begin
        p_vec = [t[1] for t in col]
        T_vec = [t[2] for t in col]
        q_vec = [t[3] for t in col]
        lev_to_z_column_pTq(p_vec, T_vec, q_vec; thermodynamics_backend)
    end, pqT; dims = ldn)
end

function lev_to_z(
    p::FT,
    T::FT,
    q::FT,
    pg::FT,
    Tg::FT,
    qg::FT;
    data,
    thermodynamics_backend =DefaultThermodynamicsBackend(),
) where {FT <: Real}
    return lev_to_z([p], [T], [q], [pg], [Tg], [qg]; data, thermodynamics_backend)
end

function lev_to_z_from_LES_output(
    p,
    T,
    q,
    pg,
    Tg,
    qg;
    data,
    thermodynamics_backend =DefaultThermodynamicsBackend(),
    assume_monotonic::Bool = false,
    flight_number::Int,
    forcing_type::AbstractForcingType,
    ground_indices = nothing,
)
    dimnames = NC.dimnames(data["T"])
    lev_dim_num = findfirst(x -> x == "lev", dimnames)
    ldn = lev_dim_num
    L = size(p, lev_dim_num)
    time_dim_num = findfirst(x -> x == "time", dimnames)
    tdn = time_dim_num
    Lt = size(T, time_dim_num)   # forcing time count; `p` is the time-invariant lev grid (singleton time)

    if !assume_monotonic
        LES_data = open_atlas_les_output(flight_number, forcing_type).data
        p_LES = LES_data["PRES"]
        z_LES = LES_data["z"]
        dimnames_LES = NC.dimnames(p_LES)
        tdn_LES = findfirst(x -> x == "time", dimnames_LES)

        summary_file = atlas_socrates_summary_file(flight_number)
        SOCRATES_summary = NC.Dataset(summary_file, "r")
        flight_ind = findfirst(vec(Array(SOCRATES_summary["flight_number"])) .== flight_number)
        initial_time = SOCRATES_summary["reference_time"][flight_ind] - Dates.Hour(12)

        t_in = Array(data["tsec"])
        bdate_data = Array(data["bdate"])
        bdate_scalar = if bdate_data isa AbstractArray
            non_missing = collect(skipmissing(vec(bdate_data)))
            isempty(non_missing) && error("bdate is missing for flight $(flight_number)")
            non_missing[1]
        else
            bdate_data
        end
        bdate_str = lpad(string(Int(round(bdate_scalar))), 6, '0')
        t_base = Dates.DateTime(bdate_str, Dates.DateFormat("yymmdd")) + Dates.Year(2000)
        t_in = t_base .+ Dates.Second.(t_in)
        t_in = (t_in .- initial_time) ./ Dates.Millisecond(1)
        t_in ./= 1000

        t_les = Array(LES_data["time"])
        t_les = (t_les .- t_les[1]) .* (24 * 3600)
        t_les = NCDatasets.nomissing(t_les)

        p_LES = NC.nomissing(Array(p_LES)) .* 100
        z_LES = NC.nomissing(Array(z_LES))

        p_LES = mapslices(
            x -> Interpolation.interpolate_1d(t_in, t_les, x, Interpolation.FastLinear1DInterpolation; bc = Interpolation.NearestBoundaryCondition()),
            p_LES;
            dims = 2,
        )

        FTLES = eltype(z_LES)
        s_z = collect(size(T))   # z spans the forcing times (`p` is the time-invariant lev grid); +1 lev for the ground
        s_z[ldn] += 1
        z = Array{FTLES}(undef, s_z...)

        increasing_p = false
        len = length(selectdim(p_LES, tdn_LES, 1))
        p_LES_t = Vector{FTLES}(undef, len)

        len = length(selectdim(p, tdn, 1))
        p_t_no_g = Vector{FTLES}(undef, len)
        T_t_no_g = Vector{Float64}(undef, len)
        q_t_no_g = Vector{Float64}(undef, len)
        p_t = Vector{FTLES}(undef, len + 1)
        T_t = Vector{Float64}(undef, len + 1)
        q_t = Vector{Float64}(undef, len + 1)

        for i_t in 1:Lt
            p_t_no_g .= vec(Array(selectdim(p, tdn, min(i_t, size(p, tdn)))))   # `p` is time-invariant (isobaric lev grid)
            T_t_no_g .= vec(Array(selectdim(T, tdn, i_t)))
            q_t_no_g .= vec(Array(selectdim(q, tdn, i_t)))
            pg_i = pg[i_t]
            Tg_i = Tg[i_t]
            qg_i = qg[i_t]

            p_LES_t .= vec(Array(selectdim(p_LES, tdn_LES, i_t)))
            p_LES_min, p_LES_max = extrema(p_LES_t)
            p_s_in = Float64(pg_i)

            if isnothing(ground_indices)
                index = searchsortedfirst(p_t_no_g, p_s_in)
            else
                index = only(vec(Array(selectdim(ground_indices, tdn, i_t))))
            end

            increasing_p = p_t_no_g[1] < p_t_no_g[end]
            if increasing_p
                sort!(p_LES_t)
                sort!(z_LES, rev = true)
            else
                reverse!(p_t_no_g)
                reverse!(T_t_no_g)
                reverse!(q_t_no_g)
                sort!(p_LES_t)
                sort!(z_LES, rev = true)
                index = length(p_t) - index + 1
            end

            n_work = len + 1
            if index > 1
                p_t[1:(index - 1)] .= p_t_no_g[1:(index - 1)]
                T_t[1:(index - 1)] .= T_t_no_g[1:(index - 1)]
                q_t[1:(index - 1)] .= q_t_no_g[1:(index - 1)]
            end
            p_t[index] = p_s_in
            T_t[index] = Float64(Tg_i)
            q_t[index] = Float64(qg_i)
            if index <= len
                p_t[(index + 1):n_work] .= p_t_no_g[index:end]
                T_t[(index + 1):n_work] .= T_t_no_g[index:end]
                q_t[(index + 1):n_work] .= q_t_no_g[index:end]
            end

            i_t_min_p = findfirst(x -> x > p_LES_min, p_t)
            i_t_max_p = findlast(x -> x < p_LES_max, p_t)

            selectdim(z, tdn, i_t)[i_t_min_p:i_t_max_p] =
                lev_to_z_from_LES_output_column(p_t[i_t_min_p:i_t_max_p], z_LES, p_LES_t)

            new_z = lev_to_z_column_pTq(p_t[i_t_max_p:end], T_t[i_t_max_p:end], q_t[i_t_max_p:end]; thermodynamics_backend)
            new_dz = new_z[2:end] .- new_z[1:(end - 1)]
            new_dz = cumsum(new_dz)
            new_z = selectdim(z, tdn, i_t)[i_t_max_p] .+ new_dz
            selectdim(z, tdn, i_t)[(i_t_max_p + 1):end] = new_z

            selectdim(z, tdn, i_t)[1:(i_t_min_p - 1)] =
                lev_to_z_column_pTq(p_t[1:i_t_min_p], T_t[1:i_t_min_p], q_t[1:i_t_min_p]; thermodynamics_backend)[1:(end - 1)] .+
                selectdim(z, tdn, i_t)[i_t_min_p]
            z_col = selectdim(z, tdn, i_t)
            z_col .-= z_col[index]
        end

        if !increasing_p
            z = reverse(z, dims = ldn)
        end
        return z
    else
        error("not implemented")
    end
end

function lev_to_z_from_LES_output_column(
    p_col::AbstractArray{<:Real},
    lesz,
    lesp;
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
)
    p_in = Float64.(p_col)
    return Interpolation.interpolate_1d(p_in, lesp, lesz, interp_method; bc = Interpolation.ErrorBoundaryCondition())
end

function z_from_data(data; thermodynamics_backend =DefaultThermodynamicsBackend())
    T = data["T"]
    q = data["q"]
    p = align_along_dimension(data["lev"], dim_num("lev", T))
    pg = data["Ps"]
    Tg = data["Tg"]
    qg = saturation_q_tot_from_pgTg.(Ref(thermodynamics_backend), pg, Tg)
    return lev_to_z(p, T, q, pg, Tg, qg; data, thermodynamics_backend)
end

