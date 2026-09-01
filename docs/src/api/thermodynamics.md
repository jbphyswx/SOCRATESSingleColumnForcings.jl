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

## Latent heats

```@docs
latent_heat_generic
latent_heat_vapor
latent_heat_sublim
```

## Constants

Each is a backend accessor, so an extension backend supplies its own parameter set's value.

```@docs
R_d
R_v
grav
cp_d
cp_v
cp_l
cp_i
molmass_ratio
p_ref
T_0
T_freeze
T_icenuc
L_v0
L_s0
e_ref
```
