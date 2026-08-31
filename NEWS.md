# NEWS

## v0.16.0

**Breaking**, and two fixes move numbers.

### Numbers change

- **Units of LES-archive variables are corrected.** `src/units.jl` now classifies every one of the
  492 variables the archive carries, recording the authority for each classification. Corrections
  against the conversions previously in use: the 34 microphysical **number**-process rates are per
  **kg**, not per cm³ (were off by 1e6); the 43 **mass**-process rates are already `kg/kg/s` in M2005,
  so no `g→kg` applies (off by 1000); `TQ` and the `QW*` flux-budget terms likewise (1000);
  `Q2`/`QTSTOR` are moisture rates, not temperature rates; `ZCB2` is a height variance, not a height;
  `PREC2` is a squared unit and takes the squared factor; `LHFOBS`/`SHFOBS` carry an **undeclared**
  `-99999` sentinel (the declared `missing_value` is `-9999`) that is now masked to `NaN`.
- **`:dTdt_rad` moves 150 s later relative to model time 0.** LES statistics are means labelled at
  their averaging-window *midpoint* — SAM stamps `day - nstat*dt/2/86400`, with `nstat` the steps per
  write — so the record opens half a period before its first label. Atlas-input fields were already
  referenced to the run start, so the two sources previously sat half a period apart inside a single
  returned `NamedTuple`; they now share an origin (`regrid_source_t_origin`).

### Fixed

- **Conservative regridding could silently use another field's mass matrix.** The cache was keyed on
  the interpolation method type alone, while the matrix also depends on the grid, `k` and the boundary
  condition — so with a per-field `new_z`, fields on different grids shared one matrix. Replaced by
  `MassMatrixCache`, keyed on all four.
- **Surface `qg` used the offset temperature instead of the SST** on flights with a positive
  `Tg_offset` (RF01, RF10): `Tg_orig` aliased an array that was then mutated in place.
- **PCHIP supported only one boundary condition** — `build_spline` threw for `Error` and rejected
  `Nearest`. It now honours all three; `PCHIPExtrapolatedInterpolant` is renamed `PCHIPInterpolant`.
- Conservative regridding now works for the **Interpolations** backend: `safe_integrate` splits at the
  knots and applies a Gauss–Legendre rule matched to the spec's polynomial degree, so it is exact.
- `coerce_to_shared_nodes` works for **all four** backends, via the new `interpolant_nodes` /
  `rebuild_interpolant` contract verbs; the `Vector`/`SVector`/`NTuple` methods collapsed into one.
- `UniformRange` is constructible on an integer-seconds axis — the package's own LES time axis.
- The docs build no longer depends on an absolute path to a worktree that does not exist.

### Changed

- `get_column_forcing` takes `z_regrid_opts::RegriddingOpts` and
  `t_regrid_opts::RegriddingOpts` in place of `interp_method` / `interp_kwargs` /
  `conservative_interp` / `conservative_interp_kwargs` / `bc` / `drop_collinear`. Conservative
  regridding is selected by *holding* a `ConservativeRegridingOpts` rather than by a separate flag, so
  the settings cannot disagree with whether it is on, and each stage carries its own boundary
  condition — which is what makes them independently settable.
- `A_cache` / `Af_cache` → a single `mass_matrix_cache::MassMatrixCache`. The factorization is now
  `lu` rather than `factorize`, making the stored type concrete; this gives up `factorize`'s choice of
  a symmetric factorization, which costs nothing here because the mass matrix divides each row by its
  own cell width and so is not symmetric unless every cell width is equal.
- `integrate_method` is `IntegrateMass()` / `InvertMass()` rather than `:integrate` / `:invert`.
- `output_interp_kwargs(Val)` → `output_z_regrid_opts(Val, opts)`.
- `les_reference_profiles` takes one `z_regrid_opts`.
- Added `q_vap_saturation_ice` to the thermodynamics contract and the Thermodynamics extension.

### Performance

Allocation work in `get_column_forcing` (baseline 18.46 MiB for 11 fields, **not yet re-measured**):
slice copies replaced by views, the shared coordinate hoisted out of the per-column loop,
`combine_air_and_ground_data`'s insert path written into a preallocated output instead of `mapslices`,
`lev_to_z` no longer building an array of tuples and its column kernel allocating nothing, and the
evaluate path preallocating its output and using an in-place kernel where the method has one.

## v0.15.1

- Fixed RF09's default vertical grid to use its native 320-level LES height coordinate
  (`25–6047.79 m`) instead of the shorter shared 320-level grid used by RF01, RF10, and RF11.
- The `atlas_les_metadata_v1` artifact now includes `RF09_grd.txt`; other flights continue to use
  the shared `192level-grd.txt` or `320level-grd.txt`.

## v0.15.0

Major source reorganization and performance-oriented interpolation storage. **Breaking** for code that called `process_case`, unqualified interpolation symbols, or paths under removed files.

### Thermodynamics

- Supports `Thermodynamics = "1"` (via the package extension).

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

- Lowered `Integrals.jl` compatibility to keep compatibility with older downstream stacks.

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

- Updated to Thermodynamics 0.11 for newer downstream integration.

## v0.9.15

- Updated surface extrapolation logic toward the ground SST / Tg state.

## v0.9.14

- Added test updates and dependency refreshes around LES-output-based forcing work.

## v0.9.13

- Added LES-output-based `z` handling support, Atlas output readers, and Tg-offset logic refinements.

## v0.1.0 - v0.9.12

- Early package development is tracked in the repository tags and GitHub-generated release notes.
- A later pass can expand these entries if you want a completely reconstructed changelog from the old tag messages.