# Forcings API

```@meta
CurrentModule = SSCF
```

## Forcing types

```@docs
AbstractForcingType
ObsForcing
ERA5Forcing
forcing_key
symbol
flight_numbers
forcing_types
grid_height
is_valid_flight_number
```

## Column forcing

```@docs
supported_forcing_variables
get_column_forcing
default_new_z
les_reference_profiles
les_reference_profiles!
```

## Regridding options

The two regrid stages are configured independently; see [Column forcings](../forcings.md) for how they
compose.

```@docs
RegriddingOpts
RegriddingOpts
ConservativeRegridingOpts
output_z_regrid_opts
MassMatrixCache
MassMatrixKey
mass_matrix!
```

## LES averaging windows

```@docs
les_output_period
les_averaging_window
```

## Surface

```@docs
get_surface_reference_state
get_surface_reference_state!
get_surface_forcing
```
