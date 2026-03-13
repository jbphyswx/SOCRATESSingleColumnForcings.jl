# SOCRATESSingleColumnForcings

| **DOI**                           | [![DOI][zenodo-img]][zenodo-latest-url]          |
|-----------------------------------|--------------------------------------------------|

[zenodo-img]: https://zenodo.org/badge/585317234.svg
[zenodo-latest-url]: https://doi.org/10.5281/zenodo.14945665


This repository includes scripts to take in SOCRATES flight data and create both 2D [lev,time] datasets and 1D [lev] columns of time-interpoloation functions to force the CliMA EDMF Turbulence-Convection single-column model.

The profiles are created off the data used to force LES simulations by Atlas, 2020 (https://doi.org/10.1029/2020MS002205), some of her scripts are included in /Rachel_Atlas_Scripts.

The output data are too large to store here and so will be stored in box. They 1D functions are interpolated to the grid that Atlas used, but functions are provided to interpolate to your own grid/time and return either raw data in time or a column of spline functions...


Data retrieval is now artifact-backed in `v0.14`.

- Call `SSCF.download_atlas_les_inputs(...)` and/or `SSCF.download_atlas_les_outputs(...)` to create/update per-flight artifacts in your Julia artifact store.
- Call `SSCF.open_atlas_les_input(...)` and `SSCF.open_atlas_les_output(...)` to read forcing/profile files from those artifacts.
- `SOCRATES_summary.nc` is also resolved from the artifact path used by SSCF.

The raw download/link definitions remain in `/Data/Atlas_LES_Profiles/`.

There are a variety of helper functions in `src/helper_functions.jl`

However, the main setup for a case comes from running `process_case()` in `src/process_case.jl`

This constructs the forcing for your chosen flight on your chosen new vertical grid and time grids by interpolating Atlas's data onto the vertical grid, and, depending on if you chose to return the initial conidtion or not, either returns that initial conditoin or converts the time dimension into splines which can be evaluated at your chosen timestep. 


The model returns forcings ` dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge`

`H_nudge` and `qt_nudge` are taken from your choice of forcings and represent relaxation of the thermal profile (liquid-ice potential temperature) and total moisture respectively

The remainder of the properties are always taken from the ERA5 based Atlas inputs.



Once this is finalized, perhaps the relevant parts can be ported to `CliMA/AtmosphericProfilesLibrary.jl` (`CliMA/TurbulenceConvection.jl` as is doesn't access this repo)
