#=
Interpolation submodule for SOCRATESSingleColumnForcings.

Public surface is reached through qualified calls: `Interpolation.foo` from package source and
`SSCF.Interpolation.foo` from extensions (no names are exported into the parent namespace).

Package-source callers use `Interpolation.<name>`; extensions extend/subtype
`SSCF.Interpolation.<name>` (e.g. `build_spline`, `AbstractInterpolationMethod`, `nnls_solve`).

Include order matters: boundary-condition types first, then the interpolants that reference them,
then the conservative helpers that use `build_spline`/`safe_integrate`.
=#
module Interpolation

using LinearAlgebra: LinearAlgebra
using Statistics: Statistics
using StaticArrays: StaticArrays

include("boundary_conditions.jl")
include("interpolants.jl")
include("conservative.jl")

end
