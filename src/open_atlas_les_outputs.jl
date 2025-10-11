#

"""
open_atlas_les_output(flight_number::Int, forcing_type::Symbol; open_files::Bool = true, include_grid::Bool = true)

opens the files downloaded in download_atlas_les_profiles.jl
"""
function open_atlas_les_output(
    flight_number::Int,
    forcing_type::Symbol;
    open_files::Bool = true,
    include_grid::Bool = true,
)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles", "Output_Data") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number, pad = 2)

    if forcing_type == :obs_data
        data_filename = joinpath(atlas_dir, RF_num * "_Obs_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    elseif forcing_type == :ERA5_data
        data_filename = joinpath(atlas_dir, RF_num * "_ERA5_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc
    else
        error("forcing_type must be :obs_data or :ERA5_data")
    end


    # check if file exists and if not, download it
    # The data is quite large so we'll try to load it from Box
    if !isfile(data_filename)
        download_atlas_les_outputs(; flight_numbers = (flight_number,), forcing_types = (forcing_type,)) # Atlas website
    end


    data = isfile(data_filename) ? (open_files ? NC.Dataset(data_filename, "r") : data_filename) : nothing

    if include_grid
        grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")
        grid_data = isfile(grid_filename) ? (open_files ? vec(readdlm(grid_filename, FT)) : grid_filename) : nothing
        return NamedTuple{(forcing_type, :grid_data)}((data, grid_data))
    else
        return NamedTuple{(forcing_type,)}((data,))
    end
end


function open_atlas_les_output(flight_number::Int; open_files::Bool = true, include_grid::Bool = true)
    FT = Float64 # idk how to pass on a type here without necessarily having to give a variable...
    atlas_dir = joinpath(dirname(@__DIR__), "Data", "Atlas_LES_Profiles", "Output_Data") # doesn't seem to work to use @__DIR__ directly as a variable
    RF_num = "RF" * string(flight_number, pad = 2)
    obs_filename = joinpath(atlas_dir, RF_num * "_Obs_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_obs-based_SAM_input.nc
    ERA5_filename = joinpath(atlas_dir, RF_num * "_ERA5_" * "SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc") # e.g. https://atmos.uw.edu/~ratlas/RF12_ERA5-based_SAM_input_mar18_2022.nc

    # check if file exists and if not, download it
    # The data is quite large so we'll try to load it from Box
    if !isfile(obs_filename)
        download_atlas_les_outputs(; flight_numbers = (flight_number,), forcing_types = (:obs_data,)) # Atlas website
    end
    if !isfile(ERA5_filename)
        download_atlas_les_outputs(; flight_numbers = (flight_number,), forcing_types = (:ERA5_data,)) # Atlas website
    end

    obs_data = isfile(obs_filename) ? (open_files ? NC.Dataset(obs_filename, "r") : obs_filename) : nothing
    ERA5_data = isfile(ERA5_filename) ? (open_files ? NC.Dataset(ERA5_filename, "r") : ERA5_filename) : nothing

    if include_grid
        grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")
        grid_data = isfile(grid_filename) ? (open_files ? vec(readdlm(grid_filename, FT)) : grid_filename) : nothing
        return (; obs_data, ERA5_data, grid_data)
    else
        return (; obs_data, ERA5_data)
    end
end
