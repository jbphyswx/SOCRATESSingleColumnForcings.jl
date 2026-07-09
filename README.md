# SOCRATESSingleColumnForcings.jl

| **Docs** | [![Documentation][docs-img]][docs-url] |
| **DOI** | [![DOI][zenodo-img]][zenodo-latest-url] |

[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/
[zenodo-img]: https://zenodo.org/badge/585317234.svg
[zenodo-latest-url]: https://doi.org/10.5281/zenodo.14945665

Julia package for building **single-column forcings** from the [SOCRATES](https://doi.org/10.1029/2019JD031915) field campaign, using the LES input/output datasets published by [Atlas (2020)](https://doi.org/10.1029/2020MS002205). The primary consumer is the CliMA EDMF single-column model ([TurbulenceConvection.jl](https://github.com/CliMA/TurbulenceConvection.jl)), but the package is usable on its own for reading Atlas data, regridding profiles, and building allocation-free time interpolants.

## What it does

Given a SOCRATES **flight number** and a **forcing source** (observation-based or ERA5-based Atlas inputs), the package:

1. Loads artifact-backed NetCDF forcing files and optional LES output files.
2. Reconstructs thermodynamic state from Atlas liquid-ice temperature fields.
3. Converts pressure levels to geometric altitude (hypsometric integration).
4. Regrids fields onto a target vertical grid (optionally mass-conservative).
5. Returns **per-level time interpolants** (or time-0 profiles) for nudging, advective tendencies, subsidence, winds, and LES radiation.

The return type is a **concretely typed `NamedTuple`**: each requested field is a `Vector` (one entry per vertical level) of built interpolants that evaluate allocation-free at runtime when storage is chosen appropriately.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jbphyswx/SOCRATESSingleColumnForcings.jl")
# or, for local development:
# Pkg.develop(path="/path/to/SOCRATESSingleColumnForcings.jl")
```

**Requirements:** Julia ≥ 1.10. Core deps: `NCDatasets`, `StaticArrays`, stdlibs.

**Recommended optional deps** (loaded as [package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))):

| Extension | Weakdep | Purpose |
|-----------|---------|---------|
| `SOCRATESSingleColumnForcingsThermodynamicsExt` | `Thermodynamics` | Accurate saturation adjustment, θ\_liq\_ice, density |
| `SOCRATESSingleColumnForcingsNonNegLeastSquaresExt` | `NonNegLeastSquares` | Positivity enforcement in conservative regridding |
| `SOCRATESSingleColumnForcingsPCHIPInterpolationExt` | `PCHIPInterpolation`, `ForwardDiff` | PCHIP splines + conservative PCHIP |
| `SOCRATESSingleColumnForcingsDierckxExt` | `Dierckx` | Cubic spline backend |
| `SOCRATESSingleColumnForcingsInterpolationsExt` | `Interpolations` | Interpolations.jl backend |

Load an extension by `using` its weak dependency in the same session, e.g. `using Thermodynamics` before calling forcing functions with a `ThermodynamicsParameters` backend.

## Quick start

```julia
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# 1. Download data (once per flight; stored as Julia artifacts)
SSCF.download_atlas_les_inputs(flight_numbers = [9])
SSCF.download_atlas_les_outputs(flight_numbers = [9])

# 2. Build column forcing on the default Atlas vertical grid
using Thermodynamics: Thermodynamics as TD
using ClimaParams: ClimaParams as CP
FT = Float32
tp = TD.Parameters.ThermodynamicsParameters(CP.create_toml_dict(FT))

forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing();
    thermodynamics_backend = tp,
)

# 3. Evaluate a nudging field at model time t = 3600 s, level k = 10
t = FT(3600)
k = 10
forcing.H_nudge[k](t)   # liquid-ice potential temperature [K]
forcing.subsidence[k](t) # subsidence [m/s]
```

Request only the fields you need:

```julia
forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing(),
    (:dTdt_hadv, :H_nudge, :subsidence);
    thermodynamics_backend = tp,
)
```

Regrid onto a custom vertical grid:

```julia
new_z = collect(0.0:100.0:4000.0)
forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing(), new_z;
    thermodynamics_backend = tp,
)
```

## Forcing variables

`supported_forcing_variables` lists every output `get_column_forcing` can produce:

| Symbol | Description | Source |
|--------|-------------|--------|
| `:dTdt_hadv` | Horizontal advective temperature tendency | Atlas input (`divT`) |
| `:dqtdt_hadv` | Horizontal advective moisture tendency | Atlas input (`divq`) |
| `:H_nudge` | Liquid-ice potential temperature profile (nudging target) | Derived from `T`, `p`, `q` |
| `:T_nudge` | Absolute temperature profile (nudging target) | `T` (+ surface `Tg`) |
| `:qt_nudge` | Total specific humidity profile (nudging target) | `q` (+ surface `qg`) |
| `:subsidence` | Large-scale subsidence | Derived from `omega`, `Ps` tendency |
| `:u_nudge`, `:v_nudge` | Horizontal wind components | Atlas input |
| `:ug_nudge`, `:vg_nudge` | Geostrophic wind components | Atlas input |
| `:dTdt_rad` | Radiative temperature tendency | LES output (`RADQR`) |

**Forcing source:** Pass `ObsForcing()` or `ERA5Forcing()` to select observation-based vs ERA5-based Atlas input files. Winds, subsidence, and geostrophic components are always ERA5-sourced regardless of `forcing_type`; nudging targets follow the chosen source.

**Flights:** `flight_numbers = (1, 9, 10, 11, 12, 13)`. RF11 is ERA5-only; RF12/RF13 use 192-level grids; others use 320 levels.

See [docs/forcings.md](docs/src/forcings.md) for the full `get_column_forcing` API, or the [Documenter site](https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/forcings/).

## Fast interpolation storage

Built time interpolants store their node/value arrays in caller-selected backings via type parameters:

```julia
# Default: StepRangeLen time axis + Vector{Float64} values
SSCF.get_column_forcing(9, SSCF.ObsForcing(); thermodynamics_backend = tp)

# Fast uniform time axis + Float32 values (recommended for production)
SSCF.get_column_forcing(
    9, SSCF.ObsForcing(),
    SSCF.supported_forcing_variables,
    Tuple{SSCF.Interpolation.UniformRange, Float32},
    Tuple{Vector, Float32};
    thermodynamics_backend = tp,
    drop_collinear = Val(true),
)
```

| Storage type | Role | When to use |
|--------------|------|-------------|
| `StepRangeLen` | Uniform time axis (stdlib range) | Default; O(1) interval search |
| `UniformRange` | Custom uniform range with precomputed `inv_step` | ~4 ns eval; avoids `StepRangeLen` twiceprecision |
| `Vector` / `SVector` | Irregular or small fixed node sets | General backing; `SVector` for isbits hot paths |
| `Constant` / `ConstantVector` | Exactly-constant fields after collinear pruning | Constant-field fast path |

All interpolation API lives under `SSCF.Interpolation.*` (qualified calls only — nothing is re-exported). See [docs/interpolation.md](docs/src/interpolation.md) or the [Documenter guide](https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/interpolation/).

## Other entry points

```julia
# Surface reference state at t = 0
SSCF.get_surface_reference_state(9, SSCF.ObsForcing(); thermodynamics_backend = tp)

# Time-dependent surface conditions as built interpolants
SSCF.get_surface_forcing(9, SSCF.ObsForcing(); thermodynamics_backend = tp)

# LES reference pressure/density profiles for TurbulenceConvection setup
SSCF.les_reference_profiles(9; forcing_type = SSCF.ObsForcing())

# Low-level I/O
SSCF.open_atlas_les_input(9, SSCF.ObsForcing())
SSCF.open_atlas_les_output(9, SSCF.ObsForcing())
SSCF.open_atlas_les_grid(9)
```

## Data and artifacts

Forcing and LES files are too large for the repository. They are fetched via `download_atlas_les_inputs` / `download_atlas_les_outputs` and registered as [Julia artifacts](https://pkgdocs.julialang.org/v1/artifacts/) (`Artifacts.toml`). Raw download URLs and Rachel Atlas's original scripts live under `Data/Atlas_LES_Profiles/`.

See [docs/data-and-artifacts.md](docs/src/data-and-artifacts.md) or the [Documenter guide](https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/data-and-artifacts/).

## Testing

```julia
using Pkg; Pkg.activate("test"); Pkg.instantiate()
# Full suite (unit + integration; integration skips flights without local data)
julia --project=test test/runtests.jl

# Unit tests only
julia --project=test test/runtests.jl unit

# Skip integration
julia --project=test test/runtests.jl no_integration

# Single unit group
julia --project=test test/runtests.jl unit_interp
```

Integration tests are data-guarded: they skip flight/forcing combinations whose artifact files are not present locally.

## Documentation

**[Documenter site (dev)](https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/)** — built from `docs/src/` on every push to `main`.

| Guide | Source |
|-------|--------|
| Getting started | [docs/src/getting-started.md](docs/src/getting-started.md) |
| Column forcings | [docs/src/forcings.md](docs/src/forcings.md) |
| Interpolation | [docs/src/interpolation.md](docs/src/interpolation.md) |
| Data and artifacts | [docs/src/data-and-artifacts.md](docs/src/data-and-artifacts.md) |
| Release notes | [NEWS.md](NEWS.md) |

Build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
# open docs/build/index.html
```

## Project layout

```
src/
  SOCRATESSingleColumnForcings.jl   # module entry, forcing types, artifact paths
  forcings.jl                       # get_column_forcing, surface conditions
  regrid.jl                         # vertical/time regridding pipeline
  field_altitude.jl                 # lev → z (hypsometric)
  ground_insertion.jl               # surface row insertion in 4-D fields
  netcdf_fields.jl                  # NetCDF read helpers
  thermodynamics.jl                 # method stubs + DefaultThermodynamicsBackend
  interpolation/                    # self-contained Interpolation submodule
  open_atlas_les_inputs.jl          # Atlas input I/O
  open_atlas_les_outputs.jl         # Atlas LES output I/O
  les_reference_profiles.jl         # p/ρ reference profiles for TC.jl
ext/                                # optional-backend extensions
Data/Atlas_LES_Profiles/            # download script + upstream links
test/                               # unit + integration tests
```

## Citation

If you use this package, please cite the SOCRATES campaign and the Atlas LES forcing dataset:

- Pfister et al. (2022), *JGR Atmospheres*, [doi:10.1029/2019JD031915](https://doi.org/10.1029/2019JD031915)
- Atlas (2020), *JAMES*, [doi:10.1029/2020MS002205](https://doi.org/10.1029/2020MS002205)

## License

See repository license file.
