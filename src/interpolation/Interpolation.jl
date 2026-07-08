#=
Interpolation submodule for SOCRATESSingleColumnForcings.

Public surface is reached through qualified calls: `Interpolation.foo` from package source and
`SSCF.Interpolation.foo` from extensions (no names are exported into the parent namespace).
=#
"""
    Interpolation

Self-contained 1D interpolation, node coercion, collinear pruning, and conservative regridding.
Call as `SOCRATESSingleColumnForcings.Interpolation.foo` — symbols are not re-exported to the
parent module.
"""
module Interpolation

using LinearAlgebra: LinearAlgebra
using Statistics: Statistics
using StaticArrays: StaticArrays

include("boundary_conditions.jl")
include("interpolants.jl")
include("conservative.jl")

end
