# Data and artifacts

SOCRATES/Atlas LES forcing files are distributed separately from this repository and registered as [Julia artifacts](https://pkgdocs.julialang.org/v1/artifacts/).

## Download entry points

```julia
using SOCRATESSingleColumnForcings
const SSCF = SOCRATESSingleColumnForcings

# Atlas SAM *input* files (forcing NetCDF + vertical grid text)
SSCF.download_atlas_les_inputs(
    flight_numbers = [1, 9, 10, 11, 12, 13],
    forcing_types = (:obs_data,),           # and/or (:ERA5_data,)
)

# Atlas LES *output* files (3D/4D LES fields incl. radiation)
SSCF.download_atlas_les_outputs(
    flight_numbers = [1, 9, 10, 11, 12, 13],
    forcing_types = (:obs_data, :ERA5_data),
)
```

Both functions populate the Julia artifact store and return `nothing`. Re-running is idempotent (skips files already present).

Implementation: `Data/Atlas_LES_Profiles/download_atlas_les_profiles.jl` (included from the main module).

## Artifact registry

`Artifacts.toml` at the repository root defines one input and one output artifact per flight:

| Artifact name | Flight |
|---------------|--------|
| `atlas_les_inputs_rf01_v1` … `atlas_les_inputs_rf13_v1` | RF01, RF09, … RF13 |
| `atlas_les_outputs_rf01_v1` … `atlas_les_outputs_rf13_v1` | RF01, RF09, … RF13 |

Resolve paths programmatically:

```julia
SSCF.atlas_les_inputs_root(9; forcing_types = (:obs_data,))
SSCF.atlas_les_outputs_root(9; forcing_types = (:obs_data,))
```

## File layout (inside an artifact)

### Input artifacts

```
Input_Data/
  RF09_grd.txt                              # vertical grid [m]
  RF09_obs-based_SAM_input.nc               # Obs forcing
  RF09_ERA5-based_SAM_input_mar18_2022.nc   # ERA5 forcing
```

### Output artifacts

```
Output_Data/
  RF09_grd.txt
  RF09_Obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc
  RF09_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc
```

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
