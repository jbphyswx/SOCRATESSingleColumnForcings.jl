# Thermodynamics API

```@meta
CurrentModule = SSCF
```

The core defines a [`DefaultThermodynamicsBackend`](@ref). Load
`Thermodynamics.jl` to activate `SOCRATESSingleColumnForcingsThermodynamicsExt` for a more fully-featured thermodynamics setup, custom pameters, and more.

## Backends

```@docs
AbstractThermodynamicsBackend
DefaultThermodynamicsBackend
```

## Thermodynamics functions 

```@docs
equilibrium_condensate
air_density
virtual_temperature
liquid_ice_pottemp
dry_pottemp
saturation_adjust_pθq
q_vap_saturation_liquid
q_vap_saturation
liquid_fraction
calc_qg_from_pgTg
saturation_q_tot_from_pgTg
```
