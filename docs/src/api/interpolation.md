# Interpolation API

```@meta
CurrentModule = SSCF.Interpolation
```

Qualified access from user code: `SSCF.Interpolation.foo`.

## Module

```@docs
Interpolation
```

## Methods and interpolants

```@docs
AbstractInterpolationMethod
FastLinear1DInterpolation
AbstractInterpolant
Fast1DLinearInterpolant
build_spline
interpolate_1d
coerce_vector
```

## Storage types

```@docs
UniformRange
Constant
ConstantVector
```

## Boundary conditions

```@docs
AbstractBoundaryCondition
ExtrapolateBoundaryCondition
ErrorBoundaryCondition
NearestBoundaryCondition
create_bc
```

## Node pruning and node sets

```@docs
drop_collinear_nodes
coerce_to_shared_nodes
interpolant_nodes
rebuild_interpolant
```

## Conservative regridding

```@docs
AbstractConservativeIntegrateMethod
IntegrateMass
InvertMass
default_conservative_interp_kwargs
conservative_regridder
conservative_mass_matrix
get_conservative_interp_kwargs
```
