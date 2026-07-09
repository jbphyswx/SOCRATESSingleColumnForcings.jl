# Data and artifacts

SOCRATES/Atlas LES forcing files are distributed separately from this repository and registered as [Julia artifacts](https://pkgdocs.julialang.org/v1/artifacts/).

## Getting the data

The Atlas files are registered as **lazy artifacts** (`Artifacts.toml`), so they download and cache automatically on first use — `open_atlas_les_input`, `get_column_forcing`, etc. fetch them via `LazyArtifacts`. No manual download step is required.

For raw retrieval into a directory of your choice (e.g. to inspect or mirror the originals from the Box / UW mirrors), use the download helpers — each takes a positional `destdir`, writes the raw NetCDF (plus the grid text, for inputs) into it, and returns `destdir`:

```julia
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

SSCF.download_atlas_les_inputs("/path/to/dir"; flight_numbers = [1, 9], forcing_types = (:obs_data, :ERA5_data))
SSCF.download_atlas_les_outputs("/path/to/dir"; flight_numbers = [1, 9])
```

Implementation: `Data/Atlas_LES_Profiles/download_atlas_les_profiles.jl` (included from the main module).

## Artifact registry

`Artifacts.toml` at the repository root defines one **shared metadata** artifact plus one input and one output artifact per `(flight, forcing)`:

| Artifact name | Contents |
|---------------|----------|
| `atlas_les_metadata_v1` | `SOCRATES_summary.nc` + level grids (`192level-grd.txt`, `320level-grd.txt`) |
| `atlas_les_inputs_rfNN_{obs,era5}_v1` | one SAM input `.nc` (per flight × forcing) |
| `atlas_les_outputs_rfNN_{obs,era5}_v1` | one LES output `.nc` (per flight × forcing) |

Flights are RF01, RF09, RF10, RF11 (ERA5 only), RF12, RF13. Resolve on-disk artifact directories programmatically (lazily downloaded on first call):

```julia
SSCF.atlas_les_inputs_root(9, SSCF.ObsForcing())    # dir holding RF09's obs input .nc
SSCF.atlas_les_outputs_root(9, SSCF.ObsForcing())
SSCF.atlas_les_metadata_root()                      # dir holding SOCRATES_summary.nc + level grids
```

## File layout (inside an artifact)

Each per-`(flight, forcing)` artifact holds a single `.nc` at its root:

```
atlas_les_inputs_rf09_obs_v1/    RF09_obs-based_SAM_input.nc
atlas_les_inputs_rf09_era5_v1/   RF09_ERA5-based_SAM_input_mar18_2022.nc
atlas_les_outputs_rf09_obs_v1/   RF09_Obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc
```

The shared metadata artifact holds the flight summary and both level grids:

```
atlas_les_metadata_v1/   SOCRATES_summary.nc   192level-grd.txt   320level-grd.txt
```

A flight's vertical grid is `$(grid_heights[flight])level-grd.txt` (RF01/09/10/11 → 320, RF12/13 → 192).

## Opening datasets

```julia
# Single forcing source
inp = SSCF.open_atlas_les_input(9, SSCF.ObsForcing())
# inp.data   — NCDataset (or path if open_files=false)
# inp.grid_data — vertical grid vector (if include_grid=true)

out = SSCF.open_atlas_les_output(9, SSCF.ObsForcing())

# Both sources at once
both_in = SSCF.open_atlas_les_input(9)   # obs_data, ERA5_data, grid_data
both_out = SSCF.open_atlas_les_output(9)

# Grid only
grid = SSCF.open_atlas_les_grid(9)
```

Keyword `open_files = false` returns file paths instead of opened datasets (useful for existence checks).

## Forcing keys vs labels

| Concept | API |
|---------|-----|
| Dispatch type | `ObsForcing()`, `ERA5Forcing()` |
| Artifact/data key | `forcing_key(ft)` → `:obs_data` or `:ERA5_data` |
| Short label | `symbol(ft)` → `:Obs` or `:ERA5` |

## Data provenance

- **SOCRATES campaign:** Pfister et al. (2022), [doi:10.1029/2019JD031915](https://doi.org/10.1029/2019JD031915)
- **Atlas LES setup & files:** Atlas (2020), [doi:10.1029/2020MS002205](https://doi.org/10.1029/2020MS002205)
- **Upstream hosting:** [UW Atlas SOCRATES LES cases](https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html)
- **Download mirrors:** Caltech Box links in `download_atlas_les_profiles.jl` (`SOCRATES_LES_inputs_Box_links`, `SOCRATES_LES_outputs_Box_links`)

Original Rachel Atlas scripts referenced in the upstream README are under `Data/Atlas_LES_Profiles/` (not required for artifact-backed workflow).

## NetCDF read conventions

All reads materialize lazy NCDatasets arrays at boundaries:

```julia
arr = Array(var)           # eager materialize
clean = NC.nomissing(arr)  # strip missing → NaN (type-matched)
```

Profile helpers in `netcdf_fields.jl`:

- `read_profile_at_time(var, z_dim, time_dim, time_index)`
- `read_profiles_over_time(var, z_dim, time_dim; time_indices=...)`

See `test/unit_ncdatasets.jl` and `test/unit_shape_contracts.jl` for contracts.

## SOCRATES summary file

`SOCRATES_summary.nc` (flight metadata, SST offsets) is resolved from the artifact path layer when needed. Flight-specific SST offsets used in surface moisture logic are also hard-coded in `get_Tg_offset` (Atlas table values).

## Manual / offline data

If artifacts cannot be downloaded (firewall, HPC without outbound HTTP):

1. Download files manually from the Box links in `download_atlas_les_profiles.jl`.
2. Place them in the artifact directory layout above.
3. Or bind a local tree with `Artifacts.bind_artifact!` (advanced; see Julia artifact docs).

Integration tests skip flights whose files are absent — unit tests do not require data.
