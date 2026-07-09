# NEWS

## v0.14.0

Major source reorganization and performance-oriented interpolation storage. **Breaking** for code that called `process_case`, unqualified interpolation symbols, or paths under removed files.

### Thermodynamics

- Supports `Thermodynamics = "0.11, 0.12, 0.13, 0.14, 0.15"` (via the package extension).

### Package layout

- Split monolithic `helper_functions.jl`, `process_case.jl`, `les_reader_helper.jl`, and `interpolating_methods.jl` into focused modules: `forcings.jl`, `regrid.jl`, `netcdf_fields.jl`, `field_altitude.jl`, `ground_insertion.jl`, `array_utils.jl`, `thermodynamics.jl`, `les_reference_profiles.jl`.
- Added self-contained `Interpolation` submodule (`src/interpolation/`). All interpolation calls are qualified (`Interpolation.foo`); nothing re-exported to the parent namespace.
- Added package extensions for optional backends: Dierckx, Interpolations, PCHIP, NonNegLeastSquares, Thermodynamics.

### Forcing API

- **`get_column_forcing`** replaces `process_case`. Returns a `NamedTuple` keyed by requested `forcing_variables` (default: `supported_forcing_variables`).
- Custom `forcing_variables` subsets with validation and lazy shared precompute.
- Added `:T_nudge` output.
- Surface helpers: `get_surface_reference_state`, `get_surface_forcing`.
- Forcing source selected by type: `ObsForcing()`, `ERA5Forcing()`.

### Interpolation

- **`UniformRange`**: custom exactly-uniform axis with precomputed `inv_step` for fast linear eval.
- **`Constant` / `ConstantVector`**: backing types for exactly-constant fields.
- Storage specs `Tuple{Backing, Eltype}` on `get_column_forcing` for type-stable, allocation-free returns.
- `drop_collinear_nodes` preserves per-array backing (`Vector`, `SVector`, uniform `AbstractRange`); errors instead of demoting ranges to `Vector`.
- Conservative regridding consolidated in `interpolation/conservative.jl`.

### I/O

- `open_atlas_les_inputs.jl` / `open_atlas_les_outputs.jl`: aligned `Val` dispatch and explicit errors on missing files.
- Artifact-backed downloads unchanged in workflow; see `docs/data-and-artifacts.md`.

### Tests

- New: `integration_forcings.jl`, `unit_regrid_source.jl`, `unit_thermodynamics.jl`, `allocations.jl`, `inferrability.jl`, `extensions.jl`, `quality.jl`.
- Removed: `integration_process_case.jl`.

### Documentation

- Rewrote `README.md`.
- Added Documenter.jl site (`docs/make.jl`, `docs/src/`, CI workflow `.github/workflows/docs.yml`).
- User guides: getting started, forcings, interpolation, data/artifacts.



## v0.13.13

- Fixed cached factorization types.

## v0.13.12

- Fixed outdated cache usage in conservative interpolation.

## v0.13.11

- Updated workflows and formatting.

## v0.13.10

- Switched boundary-condition objects to concrete types to improve `isbits` behavior.

## v0.13.9

- Fixed remaining `nomissing` typo fallout and a vertical interpolation bug.
- Improved support for true `SVector` usage.

## v0.13.8

- Fixed a `nomissing` typo.

## v0.13.7

- Reduced allocations and improved performance.
- Moved more code toward concrete types and `SVector`-friendly execution.

## v0.13.6

- Updated Atlas download links and download helpers after upstream site changes.

## v0.13.5

- Added the true surface temperature to outputs so sensible and latent heat fluxes can be computed correctly.

## v0.13.4

- Fixed factorization and positivity-enforcement bugs.
- Switched the default NNLS algorithm from `:fnnls` to `:pivot` for performance.

## v0.13.3

- Refined positivity handling so it is only applied where physically warranted.
- Added speed improvements in the interpolation path.

## v0.13.2

- Added a configurable NNLS tolerance.
- Fixed NNLS scaling behavior for very small inputs.

## v0.13.1

- Fixed a `Dierckx.jl` integration bug where boundary conditions were not respected.

## v0.13.0

- Grouped `conservative_interp_kwargs` for better portability.

## v0.12.0

- Removed the `Integrals.jl` dependency.

## v0.11.1

- Lowered `Integrals.jl` compatibility to keep compatibility with older TurbulenceConvection / CalibrateEDMF stacks.

## v0.11.0

- Added conservative regridding and enabled it by default.
- Updated workflows to Julia 1.10.8.

## v0.10.3

- Performed a maintenance release with formatting, type annotations, test updates, performance cleanup, and CI updates to Julia 1.10.5.

## v0.10.2

- Added radiation support, variable-specific target grids, smoother profile interpolations, and related file updates.

## v0.10.1

- Ensured the LES output download folder is created automatically when needed.

## v0.10.0

- Updated to Thermodynamics 0.11 for newer TurbulenceConvection integration.

## v0.9.15

- Updated surface extrapolation logic toward the ground SST / Tg state.

## v0.9.14

- Added test updates and dependency refreshes around LES-output-based forcing work.

## v0.9.13

- Added LES-output-based `z` handling support, Atlas output readers, and Tg-offset logic refinements.

## v0.1.0 - v0.9.12

- Early package development is tracked in the repository tags and GitHub-generated release notes.
- A later pass can expand these entries if you want a completely reconstructed changelog from the old tag messages.