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

## Surface

```@docs
get_SSCF_surface_reference_state!
get_surface_reference_state
get_surface_forcing
```
