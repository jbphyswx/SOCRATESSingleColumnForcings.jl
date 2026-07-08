# output-file name tag per forcing source (e.g. `RF12_Obs_SOCRATES_...nc`)
_output_file_tag(::ObsForcing) = "Obs"
_output_file_tag(::ERA5Forcing) = "ERA5"

const _OUTPUT_FILE_SUFFIX = "_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"

"""
    open_atlas_les_output(flight_number, forcing_type; open_files=true, include_grid=true)

Open the Atlas LES *output* dataset for `flight_number` and `forcing_type`, returning
`(; data, grid_data)` — `data` is the opened `NCDataset` (the file path if `open_files=false`).
With `include_grid=false`, returns `(; data)`.
"""
function _open_atlas_les_output(flight_number::Integer, forcing_type::AbstractForcingType, ::Val{open_files}, include_grid_val::Val{include_grid}) where {open_files, include_grid}
    artifact_dir = atlas_les_outputs_root(flight_number; forcing_types = (forcing_key(forcing_type),))
    atlas_dir = joinpath(artifact_dir, "Output_Data")
    RF_num = "RF" * string(flight_number, pad = 2)
    data_filename = joinpath(atlas_dir, RF_num * "_" * _output_file_tag(forcing_type) * _OUTPUT_FILE_SUFFIX)

    data = isfile(data_filename) ? (open_files ? NC.Dataset(data_filename, "r") : data_filename) : error("Missing outputfile $data_filename")
    include_grid || return (; data)

    grid_filename = joinpath(atlas_dir, RF_num * "_grd.txt")
    grid_data = isfile(grid_filename) ? (open_files ? vec(DelimitedFiles.readdlm(grid_filename, Float64)) : grid_filename) : error("Missing grid file $grid_filename")
    return (; data, grid_data)
end

@inline open_atlas_les_output(flight_number::Integer, forcing_type::AbstractForcingType; open_files::Bool = true, include_grid::Bool = true) =
    _open_atlas_les_output(flight_number, forcing_type, Val(open_files), Val(include_grid))


"""
    open_atlas_les_output(flight_number; open_files=true, include_grid=true)

Convenience multi-source load: returns `(; obs_data, ERA5_data, grid_data)`.
"""
function _open_atlas_les_output(flight_number::Integer, open_files_val::Val{open_files}, include_grid_val::Val{include_grid}) where {open_files, include_grid}
    obs_data = _open_atlas_les_output(flight_number, ObsForcing(), open_files_val, Val(false)).data
    ERA5_data = _open_atlas_les_output(flight_number, ERA5Forcing(), open_files_val, Val(false)).data
    include_grid || return (; obs_data, ERA5_data)
    artifact_dir = atlas_les_outputs_root(flight_number; forcing_types = (forcing_key(ObsForcing()),))
    grid_filename = joinpath(artifact_dir, "Output_Data", "RF" * string(flight_number, pad = 2) * "_grd.txt")
    grid_data = isfile(grid_filename) ? (open_files ? vec(DelimitedFiles.readdlm(grid_filename, Float64)) : grid_filename) : error("Missing grid file $grid_filename")
    return (; obs_data, ERA5_data, grid_data)
end


@inline open_atlas_les_output(flight_number::Integer; open_files::Bool = true, include_grid::Bool = true) =
    _open_atlas_les_output(flight_number, Val(open_files), Val(include_grid))

