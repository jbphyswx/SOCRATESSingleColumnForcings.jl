# Getting started

End-to-end workflow: install the package, fetch SOCRATES/Atlas data, build column forcings, evaluate them at model runtime.

## 1. Install

```julia
using Pkg
Pkg.add(url="https://github.com/CliMA/SOCRATESSingleColumnForcings.jl")
```

For development:

```bash
git clone https://github.com/CliMA/SOCRATESSingleColumnForcings.jl.git
cd SOCRATESSingleColumnForcings.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Julia ≥ 1.10 is required.

## 2. Download data

Atlas LES input and output files are stored as Julia artifacts (see [data-and-artifacts.md](data-and-artifacts.md)).

```julia
using SOCRATESSingleColumnForcings
const SSCF = SOCRATESSingleColumnForcings

# Input forcing files (obs-based and/or ERA5-based)
SSCF.download_atlas_les_inputs(flight_numbers = [9], forcing_types = (:obs_data,))
SSCF.download_atlas_les_inputs(flight_numbers = [9], forcing_types = (:ERA5_data,))

# LES output files (needed for :dTdt_rad)
SSCF.download_atlas_les_outputs(flight_numbers = [9], forcing_types = (:obs_data, :ERA5_data))
```

Verify a file is present without opening it:

```julia
p = SSCF.open_atlas_les_input(9, SSCF.ObsForcing(); open_files = false)
isfile(p.data)  # true after download
```

## 3. Choose a thermodynamics backend

Equilibrium-derived quantities (`H_nudge`, saturation-adjusted `T`, condensate-aware density) require a thermodynamics backend.

**Accurate path** (recommended for production / TurbulenceConvection):

```julia
using Thermodynamics
using ClimaParams
const TD = Thermodynamics
const CP = ClimaParams
FT = Float32
tp = TD.Parameters.ThermodynamicsParameters(CP.create_toml_dict(FT))
```

Loading `Thermodynamics` activates `SOCRATESSingleColumnForcingsThermodynamicsExt`.

**Fallback path** (no Thermodynamics dep; naive ideal-gas physics):

```julia
tb = SSCF.DefaultThermodynamicsBackend()
```

## 4. Build column forcing

Minimal call (all supported fields, default Atlas grid, default storage):

```julia
forcing = SSCF.get_column_forcing(9, SSCF.ObsForcing(); thermodynamics_backend = tp)
```

The return value is a `NamedTuple`. Each field is a `Vector` with one built time-interpolant per vertical level.

### Custom vertical grid

```julia
new_z = collect(0.0:100.0:4000.0)  # meters
forcing = SSCF.get_column_forcing(9, SSCF.ObsForcing(), new_z; thermodynamics_backend = tp)
```

Per-field grids are also supported via a `NamedTuple` keyed by forcing variable symbols.

### Subset of fields

Only requested fields are computed (lazy precompute):

```julia
forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing(),
    (:H_nudge, :qt_nudge, :subsidence);
    thermodynamics_backend = tp,
)
```

### Initial condition instead of time splines

```julia
ic = SSCF.get_column_forcing(
    9, SSCF.ObsForcing();
    initial_condition = true,
    thermodynamics_backend = tp,
)
# ic.H_nudge[k] is a scalar (time-0 profile value), not a callable interpolant
```

### Fast storage types

For allocation-free, inferrable evaluation (~4 ns per call on uniform axes):

```julia
using StaticArrays

forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing(),
    SSCF.supported_forcing_variables,
    Tuple{SSCF.Interpolation.UniformRange, Float32},  # time-axis backing + eltype
    Tuple{StaticArrays.SVector, Float32};             # value backing + eltype
    thermodynamics_backend = tp,
    drop_collinear = Val(true),
)
```

See [interpolation.md](interpolation.md) for storage semantics.

## 5. Evaluate at runtime

When `initial_condition = false` (default), each level entry is a callable interpolant:

```julia
t = FT(7200)   # seconds since reference time
k = 15         # vertical index on new_z

θ_liq_ice = forcing.H_nudge[k](t)
w_subs    = forcing.subsidence[k](t)
```

Boundary conditions on built splines default to `ErrorBoundaryCondition()` (error outside the stored time range). Surface time series from `get_surface_conditions` use extrapolation.

## 6. Surface conditions

Reference surface state at t = 0:

```julia
surf = SSCF.get_surface_reference_state(9, SSCF.ObsForcing(); thermodynamics_backend = tp)
# surf.pg, surf.Tg, surf.q_tot_g
```

Time-dependent surface forcings (each field is a built time interpolant):

```julia
surf = SSCF.get_surface_conditions(9, SSCF.ObsForcing(); thermodynamics_backend = tp)
# surf.pg(t), surf.Tg(t), surf.Tsfc(t), surf.qg(t), surf.qsfc(t)
```

## 7. LES reference profiles (TurbulenceConvection setup)

```julia
profiles = SSCF.les_reference_profiles(9; forcing_type = SSCF.ObsForcing())
# profiles.p_c, profiles.p_f, profiles.ρ_c, profiles.ρ_f
```

## 8. Run tests

```bash
julia --project=test test/runtests.jl unit          # fast unit tests
julia --project=test test/runtests.jl               # full suite (skips missing data)
julia --project=test test/runtests.jl no_integration
```

## Common issues

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Missing input file …` | Artifact not downloaded | `download_atlas_les_inputs(...)` |
| `unsupported forcing variable` | Typo in `forcing_variables` tuple | Use symbols from `supported_forcing_variables` |
| `UniformRange requires an exactly-uniform axis` | Time axis has Float32 jitter | Use `UniformRange` only on exact grids, or use `StepRangeLen` / `Vector` |
| Naive θ\_liq\_ice values | `DefaultThermodynamicsBackend` | Load `Thermodynamics` and pass `ThermodynamicsParameters` |
| `:dTdt_rad` errors | LES output missing | `download_atlas_les_outputs(...)` |

## Next steps

- [forcings.md](forcings.md) — full `get_column_forcing` keyword reference
- [interpolation.md](interpolation.md) — storage types, conservative regrid, collinear pruning
- [data-and-artifacts.md](data-and-artifacts.md) — artifact layout and file naming
