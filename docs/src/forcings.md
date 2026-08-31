# Column forcings

This guide covers `get_column_forcing`, surface helpers, and the forcing-variable dispatch table.

## Forcing types

Forcing source is selected by singleton type (not a symbol):

```julia
SSCF.ObsForcing()   # observation-based Atlas SAM input
SSCF.ERA5Forcing()  # ERA5-based Atlas SAM input
```

| Type | Artifact key | Atlas input file tag |
|------|--------------|----------------------|
| `ObsForcing()` | `:Obs_data` | `obs-based_SAM_input.nc` |
| `ERA5Forcing()` | `:ERA5_data` | `ERA5-based_SAM_input_mar18_2022.nc` |

**Which fields follow `forcing_type`?** Nudging targets (`H_nudge`, `T_nudge`, `qt_nudge`) and advective tendencies (`dTdt_hadv`, `dqtdt_hadv`) come from the chosen input file. Horizontal winds, geostrophic winds, and subsidence are always ERA5-sourced.

## Supported flights

```julia
SSCF.flight_numbers  # (1, 9, 10, 11, 12, 13)
```

| Flight | Vertical levels | Notes |
|--------|-----------------|-------|
| 1, 9, 10, 11 | 320 | Deep boundary-layer cases |
| 12, 13 | 192 | Shallow cloud-topped cases |
| 11 | — | ERA5 forcing only (no obs-based input in artifact set) |

Grid height helper: `SSCF.grid_height(flight_number)`.

## `get_column_forcing`

```julia
get_column_forcing(
    flight_number,
    forcing_type,
    forcing_variables = supported_forcing_variables,
    interpolant_coord_types = Tuple{StepRangeLen, Nothing},
    interpolant_value_types = Tuple{Vector, Float64},
    FT = _working_eltype_from_value_spec,
    ;
    new_z = nothing,
    initial_condition = false,
    thermodynamics_backend = DefaultThermodynamicsBackend(),
    use_LES_output_for_z = false,
    return_old_z = false,
    fail_on_missing_data = true,
    z_regrid_opts = RegriddingOpts(),
    t_regrid_opts = RegriddingOpts(),
    mass_matrix_cache = MassMatrixCache{FT}(),
)
```

### Return shape

- **`initial_condition = false` (default):** `NamedTuple{forcing_variables}` where each value is `Vector{<:AbstractInterpolant}` — one time interpolant per vertical level on `new_z`.
- **`initial_condition = true`:** Same keys, but each level holds the time-0 scalar profile value (not callable).
- **`return_old_z = true`:** Returns the source altitude field `z_old` instead of forcing (ignores `forcing_variables` for the return).

### Key keyword arguments

| Keyword | Default | Description |
|---------|---------|-------------|
| `new_z` | Atlas grid | Target vertical grid (`AbstractVector`) or per-field `NamedTuple` |
| `initial_condition` | `false` | Return time-0 profiles instead of time splines |
| `thermodynamics_backend` | `DefaultThermodynamicsBackend()` | Pass `ThermodynamicsParameters` for accurate physics |
| `z_regrid_opts` | `RegriddingOpts()` | Method, boundary condition and conservative recipe for the regrid onto `new_z` |
| `t_regrid_opts` | `RegriddingOpts()` | Method, boundary condition and collinear pruning for the returned time splines |
| `use_LES_output_for_z` | `false` | Derive altitude from LES output instead of hypsometric integration |
| `mass_matrix_cache` | `MassMatrixCache{FT}()` | Reuse conservative mass matrices across fields/calls |

### The two regrid stages

Regridding evaluates onto `new_z` and then builds time splines. The two ask different questions, so
each takes its own options object rather than sharing one bag of keywords:

| | `z_regrid_opts::RegriddingOpts` | `t_regrid_opts::RegriddingOpts` |
|---|---|---|
| stage | evaluating onto `new_z` | building the **returned** interpolants |
| `method` | `FastLinear1DInterpolation` | `FastLinear1DInterpolation` |
| `bc` answers | a target level outside the source column | a model time outside the forcing record |
| `bc` default | `ErrorBoundaryCondition()`, or `Extrapolate` when conservative | `ErrorBoundaryCondition()` |
| also carries | `conservative` (see below) | `drop_collinear` — prune collinear nodes (exact, size-reducing) |

```julia
forcing = SSCF.get_column_forcing(
    9, SSCF.ObsForcing();
    z_regrid_opts = SSCF.RegriddingOpts(; bc = SSCF.Interpolation.ErrorBoundaryCondition()),
    t_regrid_opts = SSCF.RegriddingOpts(; bc = SSCF.Interpolation.NearestBoundaryCondition()),
)
```

Conservative regridding is selected by *holding* a `ConservativeRegridingOpts` rather than by a
separate flag, so the settings cannot disagree with whether it is on:

```julia
z_regrid_opts = SSCF.RegriddingOpts(;
    conservative = SSCF.ConservativeRegridingOpts(; k = 1, f_enhancement_factor = 5),
)
```

`get_surface_forcing` takes a single `bc` (it builds time interpolants only, with no vertical stage),
defaulting to `ExtrapolateBoundaryCondition()`.

### Storage type parameters

The 4th and 5th positional arguments control built interpolant storage (see [interpolation.md](interpolation.md)):

```julia
# coord spec: Tuple{Backing, Eltype} — Nothing = passthrough
# value spec: Tuple{Backing, Eltype} — eltype sets working FT

get_column_forcing(
    fl, ft,
    (:H_nudge,),
    Tuple{SSCF.Interpolation.UniformRange, Float32},
    Tuple{StaticArrays.SVector, Float32};
    thermodynamics_backend = tp,
)
```

Working float type `FT` is inferred from the value spec's eltype (`Nothing` → `Float64`).

### Lazy computation

Shared column precompute (Atlas fields: `T`, `p`, `q`, surface values, density, altitude) runs only if at least one requested variable needs it. LES load for `:dTdt_rad` runs only when that symbol is requested.

Unknown symbols error immediately:

```julia
# ERROR: unsupported forcing variable :foo — supported: (...)
get_column_forcing(9, ObsForcing(), (:foo,))
```

## Forcing variables

Defined in `supported_forcing_variables`:

```julia
(:dTdt_hadv, :H_nudge, :T_nudge, :dqtdt_hadv, :qt_nudge,
 :subsidence, :u_nudge, :v_nudge, :ug_nudge, :vg_nudge, :dTdt_rad)
```

| Symbol | Physical quantity | Pre-regrid source | Notes |
|--------|-------------------|-------------------|-------|
| `:dTdt_hadv` | ∂T/∂t (horizontal advection) | `divT` | ERA5 |
| `:dqtdt_hadv` | ∂q\_t/∂t (horizontal advection) | `divq` | ERA5 |
| `:H_nudge` | Liquid-ice potential temperature θ\_liq\_ice | Derived from `T`, `p`, `q` | Positive-definite; enhanced conservative kwargs |
| `:T_nudge` | Absolute temperature T | `T` + `Tg` | |
| `:qt_nudge` | Total specific humidity q\_t | `q` + `qg` | Positive-definite |
| `:subsidence` | Large-scale subsidence w | `omega`, `Ps` tendency | Uses density weighting |
| `:u_nudge`, `:v_nudge` | Horizontal wind | `u`, `v` | Ground value 0 |
| `:ug_nudge`, `:vg_nudge` | Geostrophic wind | `ug`, `vg` | |
| `:dTdt_rad` | Radiative heating rate | LES `RADQR` | Requires LES output artifact; K/day → K/s |

Extending outputs: add the symbol to `supported_forcing_variables` and define `compute(::Val{:your_symbol}, base, tb)` (plus `output_source` / `output_z_regrid` if non-default).

## Surface functions

### `get_surface_reference_state`

```julia
get_surface_reference_state(flight_number, forcing_type, FT = Float64;
    thermodynamics_backend = DefaultThermodynamicsBackend())
```

Returns `(; pg, Tg, q_tot_g)` scalars at the reference timestep. Surface moisture logic follows Atlas §3 (SST vs ambient-RH branch controlled by `get_Tg_offset`).

### `get_surface_forcing`

```julia
get_surface_forcing(flight_number, forcing_type,
    interpolant_coord_types = Tuple{StepRangeLen, Nothing},
    interpolant_value_types = Tuple{Vector, Float64};
    thermodynamics_backend = DefaultThermodynamicsBackend(),
    drop_collinear = Val(false),
    bc = Interpolation.ExtrapolateBoundaryCondition())
```

Returns `(; pg, Tg, Tsfc, qg, qsfc)` — each value a built time interpolant aligned to the model clock (t = 0 at the reference time), with extrapolation boundary conditions. The `Tuple{Backing, Eltype}` storage specs (as in `get_column_forcing`) give type-stable, allocation-free interpolants; the default coordinate spec stores the time axis as a `UniformRange` for O(1) evaluation

### `get_Tg_offset`

Flight-specific SST offset (Atlas table; sign convention documented in source).

## `les_reference_profiles`

Center/face pressure and density reference profiles on given vertical grids, as `FT` vectors `(; p_c, p_f, ρ_c, ρ_f)`:

```julia
les_reference_profiles(flight_number, FT = Float64;
    forcing_type = ObsForcing(),
    new_zc = nothing, new_zf = nothing,
    z_regrid_opts = RegriddingOpts(; bc = Interpolation.NearestBoundaryCondition()))
```

Default: `new_zf = [0; grid_data]`, `new_zc` at midpoints. Requires LES output for `p`, `ρ`. The vertical regrid uses the same machinery as the forcing profiles — one `z_regrid_opts` — via `var_to_new_coord`. An in-place `les_reference_profiles!(p_c, p_f, ρ_c, ρ_f, flight_number; …)` writes into caller-provided buffers (the buffers' element type sets the output type).

## `default_new_z`

Returns the Atlas default vertical grid vector for a flight:

```julia
SSCF.default_new_z(9)  # same grid as open_atlas_les_grid
```

## Conservative regridding

When `z_regrid_opts.conservative` is set, vertical regridding uses a mass matrix solved via precomputed `A` and its factorization. `mass_matrix_cache` is keyed by everything the matrix depends on — method, boundary condition, refinement factor `k`, and the target grid.

Per-field options override the caller's `z_regrid_opts` via `output_z_regrid_opts(Val(:symbol), z_regrid_opts)`:

- `:H_nudge`, `:qt_nudge` — higher enhancement factors (preserve sharp inversions)
- `:subsidence` — gentle smoothing (factor 1)

Positivity enforcement (`enforce_positivity` in conservative kwargs) requires the `NonNegLeastSquares` extension and applies only to `:H_nudge` and `:qt_nudge`.

## Integration sketch

Typical single-column driver pattern:

1. `les_reference_profiles` → initial `p`, `ρ` on column grid
2. `get_column_forcing` → nudging / advection / subsidence interpolants
3. `get_surface_forcing` → time-varying surface BCs
4. At each timestep `t`, evaluate `forcing.H_nudge[k](t)`, etc.

Pass the same `ThermodynamicsParameters` used by the model for consistent θ\_liq\_ice and saturation adjustment.
