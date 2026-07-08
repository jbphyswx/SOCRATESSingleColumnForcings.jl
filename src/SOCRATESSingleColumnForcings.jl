"""
    SOCRATESSingleColumnForcings

Build single-column forcings from SOCRATES/Atlas LES datasets for the CliMA EDMF
single-column model. Main entry point: [`get_column_forcing`](@ref). Interpolation
primitives live in the nested [`Interpolation`](@ref) submodule (qualified calls only).
"""
module SOCRATESSingleColumnForcings

using NCDatasets: NCDatasets as NC  
using NCDatasets: NCDatasets
using Artifacts: Artifacts
using Pkg: Pkg
using LazyArtifacts: LazyArtifacts  # enables lazy (download-on-first-use) resolution of the declared Artifacts.toml artifacts
using DelimitedFiles: DelimitedFiles
using Statistics: Statistics
using LinearAlgebra: LinearAlgebra
using Dates: Dates
using Downloads: Downloads
using StaticArrays: StaticArrays


const package_root = dirname(@__DIR__)
const artifacts_toml = joinpath(package_root, "Artifacts.toml")

include("thermodynamics.jl")

resolve_nan(x::FT, val::FT = zero(FT)) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0
function resolve_nan!(x::AbstractArray{FT}, val::FT = FT(0.0)) where {FT}
    @inbounds for i in eachindex(x)
        x[i] = resolve_nan(x[i], val)
    end
end
# resolve_inf(x::FT; val::FT=FT(NaN)) where {FT} = isinf(x) ? val : x # replace inf with NaN
# resolve_not_finite(x::FT, val = FT(0.0)) where {FT} = !isfinite(x) ? FT(val) : x # replace inf and nan with 0

# --- forcing source -------------------------------------------------------------------------- #
"""
    AbstractForcingType

A SOCRATES forcing source, encoded as a singleton type so behavior is selected by dispatch
(consistent with the thermodynamics/interpolation backends). Concrete: [`ObsForcing`](@ref)
(observation-based) and [`ERA5Forcing`](@ref) (ERA5-reanalysis-based).
"""
abstract type AbstractForcingType end

"""Observation-based Atlas SAM forcing source."""
struct ObsForcing <: AbstractForcingType end

"""ERA5-based Atlas SAM forcing source."""
struct ERA5Forcing <: AbstractForcingType end

"""Supported singleton forcing sources for iteration (`ObsForcing()`, `ERA5Forcing()`)."""
const forcing_types = (ObsForcing(), ERA5Forcing())

"""Supported SOCRATES flight numbers `(1, 9, 10, 11, 12, 13)`."""
const flight_numbers = (1, 9, 10, 11, 12, 13)

"""Short label for a forcing source: `:Obs` or `:ERA5`."""
symbol(::ObsForcing) = :Obs
symbol(::ERA5Forcing) = :ERA5

"""
    forcing_key(forcing_type)

The `:obs_data` / `:ERA5_data` symbol that keys the data group and the artifact/download layer —
distinct from the `:Obs` / `:ERA5` label returned by [`symbol`](@ref).
"""
forcing_key(::ObsForcing) = :obs_data
forcing_key(::ERA5Forcing) = :ERA5_data

"""Whether `flight_number` is valid for the given forcing source."""
@inline is_valid_flight_number(::ObsForcing, flight_number::Integer) = flight_number in (1, 9, 10, 12, 13)
@inline is_valid_flight_number(::ERA5Forcing, flight_number::Integer) = flight_number in (1, 9, 10, 11, 12, 13)




# Two cases with shallow cloud-topped boundary layers, RF12 and RF13, are run on a 192-level vertical grid.
# The other four cases have clouds extending through deeper boundary layers; they are run on a 320-level vertical grid.
const grid_heights = Base.ImmutableDict(1 => 320, 9 => 320, 10 => 320, 11 => 320, 12 => 192, 13 => 192) # this might be slow idk..

"""
    grid_height(flight_number)

Return the Atlas default vertical level count for `flight_number` (320 or 192).
"""
@inline function grid_height(flight_number::T) where {T <: Integer}
    if flight_number in (1, 9, 10, 11)
        return T(320)
    elseif flight_number in (12, 13)
        return T(192)
    else
        error("Invalid flight number: $flight_number")
    end
end






# struct InterpolaionsJLMethod <: AbstractInterpolationMethod end # Has so many methods

# Interpolation submodule (self-contained: stdlib + StaticArrays only). Included before the source
# files below that call `Interpolation.*`. Nothing is exported into this parent namespace; package
# source reaches the API through qualified calls (`Interpolation.foo`).
include("interpolation/Interpolation.jl")

include("../Data/Atlas_LES_Profiles/download_atlas_les_profiles.jl") # in Data/, not src/ (static path so the module stays statically analyzable)
include("metadata.jl")
include("open_atlas_les_inputs.jl")
include("open_atlas_les_outputs.jl")
include("array_utils.jl")
include("netcdf_fields.jl")
include("ground_insertion.jl")
include("field_altitude.jl")
include("regrid.jl")
include("forcings.jl")
include("les_reference_profiles.jl")

end
