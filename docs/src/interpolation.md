# Interpolation submodule

`SOCRATESSingleColumnForcings.Interpolation` is a self-contained submodule for 1D interpolation, node coercion, collinear pruning, and conservative regridding. It has no dependency on `NCDatasets` or `Thermodynamics`.

**Calling convention:** package source and user code call `SSCF.Interpolation.foo`. Extensions call `SSCF.Interpolation.foo` when extending. Nothing is exported into the parent namespace.

```
src/interpolation/
  Interpolation.jl       # module wrapper
  boundary_conditions.jl # Error / Extrapolate / Nearest BC types
  interpolants.jl        # Fast1DLinearInterpolant, storage types, build/eval
  conservative.jl        # mass matrices, NNLS, conservative_regridder
```

## Core types

### Methods

| Singleton | Description |
|-----------|-------------|
| `FastLinear1DInterpolation` | Default piecewise-linear method (fast path) |

Optional backends (via extensions): Dierckx cubic splines, PCHIP, Interpolations.jl.

### Interpolant

```julia
itp = SSCF.Interpolation.build_spline(
    SSCF.Interpolation.FastLinear1DInterpolation,
    xp, fp;
    bc = SSCF.Interpolation.ErrorBoundaryCondition(),
    drop_collinear = true,
)
val = itp(x)
```

`Fast1DLinearInterpolant` stores `(xp, fp, bc)` where `xp` and `fp` are `AbstractVector`s (or custom backings below).

## Storage coercion

`coerce_vector(spec, v)` converts raw node/value arrays into the backing requested by a type spec `Tuple{Backing, Eltype}`:

```julia
spec = Tuple{SSCF.Interpolation.UniformRange, Float32}
xp = SSCF.Interpolation.coerce_vector(spec, t_vec)
```

| `Backing` | Behavior |
|-----------|----------|
| `Nothing` | Keep input backing; optionally coerce eltype |
| `Vector` | `Vector{Eltype}` |
| `SVector` | `SVector{N, Eltype}` (length from data) |
| `AbstractRange` / `StepRangeLen` | Exactly-uniform axis; errors if non-uniform |
| `UniformRange` | Custom uniform range with `inv_step` precomputed |
| `Constant` | Single-node constant (all values equal) |
| `ConstantVector` | Length-N constant vector (all values equal) |

`Nothing` eltype = passthrough (keep array's eltype). Explicit eltype coerces via `_resolve_eltype`.

These specs are passed to `get_column_forcing` as the 4th/5th type arguments so the returned `NamedTuple` of interpolants is concretely typed.

## UniformRange

Custom uniform range optimized for hot-path linear interpolation:

```julia
r = SSCF.Interpolation.UniformRange(0.0f0, 300.0f0, 33)
# stores start, step, length, inv_step
```

- Coercion requires **exact** uniform spacing (no tolerance by default).
- `_eval_linear(::UniformRange, …)` uses multiply + `inv_step` instead of `searchsortedlast`.
- Benchmarks on this codebase: ~4 ns eval with `SVector` values + `drop_collinear=Val(true)` vs ~8 ns for default `StepRangeLen` + full `Vector`.

Use when the time axis is nominally uniform and exactly representable (or reconstructed to exact steps upstream).

## Constant / ConstantVector

For fields that are exactly constant after collinear pruning:

```julia
# Constant: length-1 AbstractVector (single value)
# ConstantVector: length-N but all entries equal
```

Callable fast paths return `fp.value` without interval search. Intended for use as `fp` backing when `drop_collinear` collapses a flat field.

## Boundary conditions

| Type | Out-of-range behavior |
|------|----------------------|
| `ErrorBoundaryCondition()` | Throw (default for column forcing) |
| `ExtrapolateBoundaryCondition()` | Linear extrapolation |
| `NearestBoundaryCondition()` | Nearest endpoint value |

Create via `create_bc(:error)`, `create_bc(:extrapolate)`, `create_bc(:nearest)`.

## drop_collinear_nodes

Removes interior nodes that lie exactly on the line through their neighbors (redundant for piecewise-linear interpolation):

```julia
xp2, fp2 = SSCF.Interpolation.drop_collinear_nodes(xp, fp)
```

**Backing-preserving:** `x` and `y` are pruned independently — `Vector` stays `Vector`, `SVector{N}` → `SVector{keep}`, `AbstractRange` → sub-range if kept nodes remain uniform (otherwise errors; no silent demotion to `Vector`).

Also defined on built interpolants:

```julia
itp2 = SSCF.Interpolation.drop_collinear_nodes(itp)
```

Pass `drop_collinear = Val(true)` to `get_column_forcing` to prune built time splines. Safe with `Vector` backing; use with `SVector`/`UniformRange` when the pruned length is predictable.

## Conservative regridding

In `conservative.jl`:

- `conservative_regridder` — mass-conserving vertical remap
- `get_conservative_A` — build/cache mass matrix
- `default_conservative_interp_kwargs` — enhancement factors, positivity flag
- `nnls_solve` — stub in core; implemented in `NonNegLeastSquares` extension

Conservative regridding is opt-in via `conservative_interp = true` in `get_column_forcing`. Density (`ρ`) is the default interpolation weight for column fields.

## Building and evaluating manually

```julia

xp = SSCF.Interpolation.coerce_vector(Tuple{SSCF.Interpolation.UniformRange, Float32}, collect(0.0f0:300.0f0:9600.0f0))
fp = SSCF.Interpolation.coerce_vector(Tuple{Vector, Float32}, rand(Float32, length(xp)))

itp = SSCF.Interpolation.Fast1DLinearInterpolant(xp, fp; bc = SSCF.Interpolation.ExtrapolateBoundaryCondition())
itp(450.0f0)

# Vectorized
SSCF.Interpolation.interpolate_1d(query_times, xp, fp, SSCF.Interpolation.FastLinear1DInterpolation)
```

## Performance notes

| Choice | Effect |
|--------|--------|
| `UniformRange` time axis | O(1) index via multiply; avoids `StepRangeLen` twiceprecision |
| `SVector` values (small N) | Isbits, fixed-size gather at eval |
| `drop_collinear = Val(true)` | Shorter node sets for flat/near-flat fields |
| `Constant` / `ConstantVector` fp | Constant-field eval fast path |
| `conservative_interp = false` | Skip mass-matrix setup (default) |

Type stability: pass concrete type specs (`Tuple{UniformRange, Float32}`, not `Tuple{Any, Float32}`) so `get_column_forcing` returns a fully concrete `NamedTuple`.

## Extensions

Extensions add methods to `SSCF.Interpolation.build_spline`, `interpolate_1d`, etc.:

| Extension file | Adds |
|----------------|------|
| `SOCRATESSingleColumnForcingsDierckxExt.jl` | `DierckxSpline1DInterpolationMethod` |
| `SOCRATESSingleColumnForcingsPCHIPInterpolationExt.jl` | PCHIP methods + conservative PCHIP |
| `SOCRATESSingleColumnForcingsInterpolationsExt.jl` | Interpolations.jl backend |
| `SOCRATESSingleColumnForcingsNonNegLeastSquaresExt.jl` | `nnls_solve` implementation |

Load the weak dependency (`using Dierckx`, etc.) to activate.

## API discipline

The submodule intentionally uses:

- `using Module: Module` (no selective imports, no `import`)
- No `export` inside the submodule
- Qualified calls only (`Interpolation.foo` / `SSCF.Interpolation.foo`)
- No parent-level forwarding aliases

Do not call unqualified `build_spline` from package source — always qualify.
