# Units and archive registry API

```@meta
CurrentModule = SSCF
```

Every variable the Atlas archive carries is declared with the units the file holds, the units it is
converted to, and the transform between them. The file's own `units` attribute is not sufficient on
its own — it is variously absent, wrong, self-referential, or shared between quantities of different
kinds — so each spec records the authority that fixed its classification.

## Reading a variable

```@docs
read_atlas_variable
convert_to_SI_units
declared_units
```

## Specs

```@docs
AtlasVarSpec
atlas_var_specs
atlas_input_var_specs
atlas_output_var_specs
add_atlas_derived_specs!
atlas_processing_order
ATLAS_COORDINATES
```

## Conversion primitives

Mixing ratio to specific content is not a scalar factor: it depends on the local total water, and its
moments take powers of `(1 - q_t)`.

```@docs
wt_to_qt
wt_to_qt!
wc_to_qc
qc_to_wc
wc2_to_qc2
wm2_to_kgm2s
K_to_Jkg
permicron_to_perm
mask_undeclared_sentinel
ATLAS_UNDECLARED_SENTINEL
ComposeFirst
```

## Archive variable groups

Names as the archive spells them, grouped by the classification they share.

```@docs
ATLAS_TEMPERATURES
ATLAS_STATIC_ENERGIES
ATLAS_WATER_MASSES
ATLAS_WATER_PATHS
ATLAS_WATER_MOMENT2
ATLAS_WATER_MASS_RATES
ATLAS_MASS_PROCESS_RATES
ATLAS_NUMBER_PROCESS_RATES
ATLAS_SLOPE_PARAMETERS
ATLAS_EFFECTIVE_RADII
ATLAS_MISLABELLED_KGKG
ATLAS_MOISTURE_RATES_IN_K
ATLAS_RADQR_KDAY
ATLAS_EMPTY_UNIT_FRACTIONS
ATLAS_ENERGY_FLUXES
ATLAS_OBSERVED_ENERGY_FLUXES
ATLAS_LIQUID_WATER_FLUXES
ATLAS_ICE_WATER_FLUXES
ATLAS_STATIC_ENERGY_MASS_FLUXES
ATLAS_TEMP_FLUX_BUDGET
ATLAS_WATER_FLUX_BUDGET
ATLAS_WATER_FLUX_BUDGET_K
ATLAS_WATER_VARIANCE_BUDGET
```

## Archive diagnostics

Quantities computed from the archive rather than read from it. They take SI arguments and return SI
results, and the constants are the ones the archive was written with.

```@docs
sedimentation_velocity
ice_distribution_slope
ice_mean_radius
ice_process_threshold_weight
ice_theory_timescale
water_vapor_diffusivity_in_air
cooper_ice_nucleus_concentration
```

```@docs
ATLAS_Q_SMALL
ATLAS_CP_D
ATLAS_LH_V0
ATLAS_LH_S0
ATLAS_SED_W_CAP
ATLAS_ICE_APPARENT_DENSITY
ATLAS_ICE_THRESHOLD_RADIUS
ATLAS_PARTICLE_MIN_RADIUS
```
