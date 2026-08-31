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
q_vap_saturation_liq
q_vap_saturation_ice
q_vap_saturation
q_vap_saturation_from_pressure
saturation_vapor_pressure
saturation_vapor_pressure_liq
saturation_vapor_pressure_ice
liquid_fraction
saturation_specific_humidity_from_pT
saturation_mixing_ratio_from_pT
```

## Phases

```@docs
AbstractPhase
Vapor
Liquid
Ice
```
