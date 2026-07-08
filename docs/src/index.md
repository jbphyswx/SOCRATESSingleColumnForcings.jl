# SOCRATESSingleColumnForcings.jl

SOCRATESSingleColumnForcings.jl builds **single-column forcings** from the [SOCRATES](https://doi.org/10.1029/2019JD031915) field campaign using the LES input/output datasets published by [Atlas (2020)](https://doi.org/10.1029/2020MS002205). The primary consumer is the CliMA EDMF single-column model ([TurbulenceConvection.jl](https://github.com/CliMA/TurbulenceConvection.jl)).

## What it does

Given a SOCRATES flight number and a forcing source (`ObsForcing()` or `ERA5Forcing()`), the package:

1. Loads artifact-backed Atlas NetCDF files.
2. Reconstructs thermodynamic state from Atlas liquid-ice temperature fields.
3. Converts pressure levels to geometric altitude.
4. Regrids onto a target vertical grid (optionally mass-conservative).
5. Returns a **concretely typed `NamedTuple`** of per-level time interpolants (or time-0 profiles).

## Quick example

```julia
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

SSCF.download_atlas_les_inputs(flight_numbers = [9])
SSCF.download_atlas_les_outputs(flight_numbers = [9])

using Thermodynamics, ClimaParams
tp = Thermodynamics.Parameters.ThermodynamicsParameters(
    ClimaParams.create_toml_dict(Float32),
)

forcing = SSCF.get_column_forcing(9, SSCF.ObsForcing(); thermodynamics_backend = tp)
forcing.H_nudge[10](3600.0f0)  # θ_liq_ice at level 10, t = 1 h
```

## User guides

| Guide | Description |
|-------|-------------|
| [Getting started](getting-started.md) | Install, download data, first forcing |
| [Column forcings](forcings.md) | `get_column_forcing`, surface functions |
| [Interpolation](interpolation.md) | Storage types, `UniformRange`, performance |
| [Data and artifacts](data-and-artifacts.md) | Downloads, file layout |

## Package highlights

- **`get_column_forcing`** — main entry point; custom `forcing_variables` subsets with lazy precompute.
- **`Interpolation` submodule** — qualified API for fast piecewise-linear splines, conservative regrid, optional backends via extensions.
- **Artifact-backed I/O** — `download_atlas_les_inputs` / `download_atlas_les_outputs` populate the Julia artifact store.
- **Type-stable returns** — caller-selected storage specs (`UniformRange`, `SVector`, `Constant`, …) for allocation-free eval.

See [References](references.md) for citations. Release notes live in [NEWS.md](https://github.com/jbphyswx/SOCRATESSingleColumnForcings.jl/blob/main/NEWS.md) on GitHub.

## API reference

Docstrings are rendered under **API** in the sidebar (`get_column_forcing`, `Interpolation.build_spline`, I/O helpers, thermodynamics seam, etc.).
