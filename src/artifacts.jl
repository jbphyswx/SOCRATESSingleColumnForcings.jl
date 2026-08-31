#=
Runtime data resolution: every Atlas input/output/metadata file the package reads is a lazy artifact
declared in `Artifacts.toml`, resolved on first use by `_artifact_root`. Nothing here writes to the
package's `Artifacts.toml`, so this works in a read-only depot and offline once cached.

Fetching the raw upstream files (to rebuild the artifact tarballs) lives in `raw_data_sources.jl`.
=#

# Artifact names: one per (flight, forcing) for the input & output forcing `.nc` files, plus a single
# shared metadata artifact holding `SOCRATES_summary.nc`, the two shared level grids, and RF09's
# distinct 320-level grid.
_rf_num(flight::Integer) = "RF" * string(flight, pad = 2)
_forcing_tag(::ObsForcing) = "obs"
_forcing_tag(::ERA5Forcing) = "era5"
_inputs_artifact_name(flight::Integer, ft::AbstractForcingType) =
    "atlas_les_inputs_" * lowercase(_rf_num(flight)) * "_" * _forcing_tag(ft) * "_v1"
_outputs_artifact_name(flight::Integer, ft::AbstractForcingType) =
    "atlas_les_outputs_" * lowercase(_rf_num(flight)) * "_" * _forcing_tag(ft) * "_v1"
const _metadata_artifact_name = "atlas_les_metadata_v1"

# Resolve a declared (lazy) artifact by name to its on-disk path, downloading on first use from the
# `[[<name>.download]]` mirror(s) in `Artifacts.toml`. `LazyArtifacts` (loaded by the module) enables
# the on-demand download.
function _artifact_root(name::AbstractString)
    Pkg.Artifacts.ensure_artifact_installed(name, artifacts_toml)
    return Pkg.Artifacts.artifact_path(Pkg.Artifacts.artifact_hash(name, artifacts_toml))
end

"""
    atlas_les_inputs_root(flight, forcing_type)

On-disk directory of the input-forcing artifact for `(flight, forcing_type)` — one `.nc` at its root
(`RFNN_<obs|ERA5>…_SAM_input…nc`). Lazily downloaded on first use (see `Artifacts.toml`).
"""
atlas_les_inputs_root(flight::Integer, forcing_type::AbstractForcingType) =
    _artifact_root(_inputs_artifact_name(flight, forcing_type))

"""
    atlas_les_outputs_root(flight, forcing_type)

On-disk directory of the output-forcing (LES) artifact for `(flight, forcing_type)` — one `.nc` at its
root. Lazily downloaded on first use.
"""
atlas_les_outputs_root(flight::Integer, forcing_type::AbstractForcingType) =
    _artifact_root(_outputs_artifact_name(flight, forcing_type))

"""
    atlas_les_metadata_root()

On-disk directory of the shared metadata artifact: `SOCRATES_summary.nc`, the shared level grids
(`192level-grd.txt`, `320level-grd.txt`), and RF09's distinct `RF09_grd.txt`.
"""
atlas_les_metadata_root() = _artifact_root(_metadata_artifact_name)

"Path to `SOCRATES_summary.nc` (in the shared metadata artifact). `flight_number` is accepted for call-site compatibility; the summary is flight-independent."
atlas_socrates_summary_file(flight_number::Integer) = joinpath(atlas_les_metadata_root(), "SOCRATES_summary.nc")

"""
    _atlas_grid_filename(flight)

Metadata filename for `flight`'s LES grid. RF09 used a distinct stretched 320-level grid; all other
flights use one of the two shared grids selected by [`grid_height`](@ref).
"""
_atlas_grid_filename(flight::Integer) =
    flight == 9 ? "RF09_grd.txt" : string(grid_heights[flight]) * "level-grd.txt"

"Path to the LES grid file for `flight` in the shared metadata artifact."
atlas_grid_file(flight::Integer) = joinpath(atlas_les_metadata_root(), _atlas_grid_filename(flight))

function _load_grid(flight_number::Integer, ::Val{open_files}, ::Type{FT} = Float64) where {open_files, FT}
    grid_filename = atlas_grid_file(flight_number)
    isfile(grid_filename) || error("Missing grid file $grid_filename")
    return open_files ? vec(DelimitedFiles.readdlm(grid_filename, FT)) : grid_filename
end
