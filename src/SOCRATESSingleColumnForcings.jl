module SOCRATESSingleColumnForcings

import Thermodynamics as TD
const TDP = TD.Parameters
import NCDatasets as NC
using DelimitedFiles
using Statistics
using Dierckx
using LinearAlgebra: LinearAlgebra, factorize
using Dates

LinearAlgebra.BLAS.set_num_threads(1) # if you're on HPC this is essential to A\b not slowing down by 5 orders of magnitude from 1ms to 100s 

resolve_nan(x::FT, val = FT(0.0)) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0
# resolve_inf(x::FT; val::FT=FT(NaN)) where {FT} = isinf(x) ? val : x # replace inf with NaN
# resolve_not_finite(x::FT, val = FT(0.0)) where {FT} = !isfinite(x) ? FT(val) : x # replace inf and nan with 0

# include our files
const FT = Float64
const flight_numbers = [1, 9, 10, 11, 12, 13]
const forcing_types = [:obs_data, :ERA5_data] # maybe change these to [:obs,:ERA5] later? would need to mirror in Cases.jl in TC.jl
const TDPS = TD.Parameters.ThermodynamicsParameters
const TDTS = TD.ThermodynamicState

# Two cases with shallow cloud-topped boundary layers, RF12 and RF13, are run on a 192-level vertical grid.
# The other four cases have clouds extending through deeper boundary layers; they are run on a 320-level vertical grid.
const grid_heights = Dict(1 => 320, 9 => 320, 10 => 320, 11 => 320, 12 => 192, 13 => 192)


const default_conservative_interp_kwargs = (;
    preserve_monotonicity = true,
    enforce_positivity = false,
    nnls_alg = :fnnls,
    nnls_tol = FT(1e-8), # default for package
    enforce_conservation = true,
    integrate_method = :invert,
)

const DCIKT = typeof(default_conservative_interp_kwargs)
const DCIKDT = Dict{Symbol, Union{Bool, Symbol, FT}} # Dict for conservative interpolation kwargs
const default_conservative_interp_kwargs_dict = DCIKDT(pairs(default_conservative_interp_kwargs))

get_conservative_interp_kwargs(::Nothing) = default_conservative_interp_kwargs

get_conservative_interp_kwargs(conservative_interp_kwargs::NamedTuple) = NamedTuple((
    key => get(conservative_interp_kwargs, key, default_conservative_interp_kwargs[key]) for
    key in keys(default_conservative_interp_kwargs)
)) # convert arbitrary NamedTuple or Dict to DCIKT, keeping only relevant keys and defaulting to the new tuple but falling back to the old one if not found

get_conservative_interp_kwargs(conservative_interp_kwargs::Dict) = get_conservative_interp_kwargs(
    NamedTuple((Symbol(k) => isa(v, String) ? Symbol(v) : v for (k, v) in conservative_interp_kwargs)),
) # convert any strings to symbols and then passed to the NamedTuple constructor


get_conservative_interp_kwargs(conservative_interp_kwargs::DCIKT) =
    merge(default_conservative_interp_kwargs, conservative_interp_kwargs) # merge between two DCIKTs (should maybe be faster than the iterative method above)

get_conservative_interp_kwargs(; kwargs...) =
    get_conservative_interp_kwargs(merge(default_conservative_interp_kwargs, kwargs)) # merge kwargs (pairs iterator) into the default DCIKT then strip it down to the relevant keys / reorder if we added any

include(joinpath("..", "Data", "Atlas_LES_Profiles", "download_atlas_les_profiles.jl")) # not in src so don't include automatically
include("open_atlas_les_inputs.jl")
include("open_atlas_les_outputs.jl")
include("helper_functions.jl")
include("process_case.jl")
# include("interpolating_methods.jl") # currently only helper_functions.jl reads these
# include("les_reader_helper.jl") # currently only helper_functions.jl reads these
include("get_LES_reference_profiles.jl")

end
