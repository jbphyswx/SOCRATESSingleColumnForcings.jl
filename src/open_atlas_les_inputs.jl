# input-file name tag per forcing source (e.g. `RF12_obs-based_SAM_input.nc`)
_input_file_tag(::ObsForcing) = "obs-based_SAM_input"
_input_file_tag(::ERA5Forcing) = "ERA5-based_SAM_input_mar18_2022"

function _open_atlas_les_input(flight_number::Integer, forcing_type::AbstractForcingType, ::Val{open_files}, ::Val{include_grid}) where {open_files, include_grid}
    data_dir = atlas_les_inputs_root(flight_number, forcing_type)   # per-(flight,forcing) artifact; .nc at its root
    RF_num = "RF" * string(flight_number, pad = 2)
    data_filename = joinpath(data_dir, RF_num * "_" * _input_file_tag(forcing_type) * ".nc")

    data = isfile(data_filename) ? (open_files ? NC.Dataset(data_filename, "r") : data_filename) : error("Missing input file $data_filename")
    include_grid || return (; data)

    grid_filename = atlas_grid_file(flight_number)   # shared metadata artifact, selected via grid_heights
    grid_data = isfile(grid_filename) ? (open_files ? vec(DelimitedFiles.readdlm(grid_filename, Float64)) : grid_filename) : error("Missing grid file $grid_filename")
    return (; data, grid_data)
end

"""
    open_atlas_les_input(flight_number, forcing_type; open_files=true, include_grid=true)

Open the Atlas LES *input* forcing dataset for `flight_number` and `forcing_type`, returning
`(; data, grid_data)` — `data` is the opened `NCDataset` (the file path if `open_files=false`).
With `include_grid=false`, returns `(; data)`.
"""
function open_atlas_les_input(flight_number::Integer, forcing_type::AbstractForcingType; open_files::Bool = true, include_grid::Bool = true)
    return _open_atlas_les_input(flight_number, forcing_type, Val(open_files), Val(include_grid))
end

function _open_atlas_les_input(flight_number::Integer, open_files_val::Val{open_files}, include_grid_val::Val{include_grid}) where {open_files, include_grid}
    obs_data = _open_atlas_les_input(flight_number, ObsForcing(), open_files_val, Val(false)).data
    ERA5_data = _open_atlas_les_input(flight_number, ERA5Forcing(), open_files_val, Val(false)).data
    include_grid || return (; obs_data, ERA5_data)
    grid_filename = atlas_grid_file(flight_number)   # shared metadata artifact
    grid_data = isfile(grid_filename) ? (open_files ? vec(DelimitedFiles.readdlm(grid_filename, Float64)) : grid_filename) : error("Missing grid file $grid_filename")
    return (; obs_data, ERA5_data, grid_data)
end

"""
    open_atlas_les_input(flight_number; open_files=true, include_grid=true)

Convenience multi-source load: returns `(; obs_data, ERA5_data, grid_data)`.
"""
function open_atlas_les_input(flight_number::Integer; open_files::Bool = true, include_grid::Bool = true)
    return _open_atlas_les_input(flight_number, Val(open_files), Val(include_grid))
end

function _open_atlas_les_grid(flight_number::Integer, ::Val{open_files}) where {open_files}
    grid_filename = atlas_grid_file(flight_number)   # shared metadata artifact, selected via grid_heights
    grid_data = isfile(grid_filename) ? (open_files ? vec(DelimitedFiles.readdlm(grid_filename, Float64)) : grid_filename) : error("Missing grid file $grid_filename")
    return (; grid_data)
end

"""
    open_atlas_les_grid(flight_number; open_files=true)

Return `(; grid_data)`, the vertical grid vector for `flight_number` (from the obs input artifact).
"""
function open_atlas_les_grid(flight_number::Integer; open_files::Bool = true)
    return _open_atlas_les_grid(flight_number, Val(open_files))
end
