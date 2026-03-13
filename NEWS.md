# NEWS

## v0.14.0

This release is the breaking compatibility port for the modern SOCRATES / ClimaAtmos calibration stack.

- Ported SSCF to `Thermodynamics = "1"`, `ForwardDiff = "1"`, and `NCDatasets = "0.14"`.
- Replaced `TD.PhaseEquil_pTq` / `TD.PhaseEquil_pθq` and the `TD.ThermodynamicState` dispatch
  type with thin internal wrappers (`phase_equil_pTq`, `phase_equil_pθq`) that return plain
  named tuples, decoupling the package from Thermodynamics internals.
- Added `air_pressure_compat`, `air_density_compat`, and `virtual_temperature_compat` helpers
  so that downstream computations on thermodynamic-state named tuples read naturally.
- Hardened all NetCDF read paths against NCDatasets v0.13+ lazy `CFVariable` /
  `DiskArrays.BroadcastDiskArray` behavior:
  - Every read boundary now calls `Array(var)` explicitly (documented NCDatasets API).
  - Added `_materialize(x)` as the single boundary guard before interpolation; uses
    `NCDatasets.nomissing(arr, FT(NaN))` with a **type-matched NaN** (via
    `Base.nonmissingtype`) so `Float32` arrays stay `Float32` rather than widening to `Float64`.
  - Removed the thin wrappers `read_all_dims` and `read_vector`; all call sites now use
    `Array(var)` / `vec(Array(var))` directly, eliminating an indirection layer that
    obscured lazy-vs-eager semantics.
- Added `read_profile_at_time` and `read_profiles_over_time` as explicit, contract-checked
  helpers that produce concrete `Vector` / `Matrix` results with shape assertions.
- Cleaned up module-level imports: no longer re-exports `nomissing`, `readdlm`, `mean`,
  `Spline1D`, `factorize`, `SVector`, or `download` into the package namespace.
- Added `package_root` and `artifacts_toml` as module-level constants replacing ad-hoc
  path constructions in download helpers.
- Added focused regression tests:
  - `unit_ncdatasets_compat.jl` — `Array(var)` shape preservation across NCDatasets versions
  - `unit_shape_contracts.jl` — profile / `(z, time)` extraction contract assertions
  - `unit_interp_methods.jl` — interpolation method correctness, including lazy/missing inputs
  - `integration_process_case.jl` — end-to-end process-case smoke test (opt-in slow path)
- Split tests into fast default unit coverage (pass `test_args=["unit","no_integration"]`)
  and an opt-in slow integration matrix.

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