# Thermodynamics API

```@meta
CurrentModule = SSCF
```

The core defines a thermo-free seam with [`DefaultThermodynamicsBackend`](@ref). Load
`Thermodynamics.jl` to activate `SOCRATESSingleColumnForcingsThermodynamicsExt` for accurate physics.

## Backends

```@docs
AbstractThermodynamicsBackend
DefaultThermodynamicsBackend
```

## Seam functions (naive backend)

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
