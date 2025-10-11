module SOCRATESSingleColumnForcings

import Thermodynamics as TD
const TDP = TD.Parameters
import NCDatasets as NC
using DelimitedFiles: DelimitedFiles, readdlm
using Statistics: Statistics, mean
using Dierckx: Dierckx, Spline1D
using LinearAlgebra: LinearAlgebra, factorize
using Dates: Dates
using Downloads: download


resolve_nan(x::FT, val::FT = FT(0.0)) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0
function resolve_nan!(x::AbstractArray{FT}, val::FT = FT(0.0)) where {FT}
    @inbounds for i in eachindex(x)
        x[i] = resolve_nan(x[i], val)
    end
end
# resolve_inf(x::FT; val::FT=FT(NaN)) where {FT} = isinf(x) ? val : x # replace inf with NaN
# resolve_not_finite(x::FT, val = FT(0.0)) where {FT} = !isfinite(x) ? FT(val) : x # replace inf and nan with 0

# include our files
const FT = Float64
const flight_numbers = (1, 9, 10, 11, 12, 13)
const forcing_types = (:obs_data, :ERA5_data) # maybe change these to [:obs,:ERA5] later? would need to mirror in Cases.jl in TC.jl
const TDPS = TD.Parameters.ThermodynamicsParameters
const TDTS = TD.ThermodynamicState

# Two cases with shallow cloud-topped boundary layers, RF12 and RF13, are run on a 192-level vertical grid.
# The other four cases have clouds extending through deeper boundary layers; they are run on a 320-level vertical grid.
const grid_heights = Dict(1 => 320, 9 => 320, 10 => 320, 11 => 320, 12 => 192, 13 => 192)


#=
For the nnls_alg,
[https://www.sciencedirect.com/science/article/pii/S1877050917307858] explains the active and passive sets well.
[https://conservancy.umn.edu/server/api/core/bitstreams/ef5b0697-b7ea-4e5b-9d97-92b29a2f468e/content] chapter 2 points out that the pivot methods tend to be computationally optimal for NNLS (Kim and Park is in NonNegLeastSquares.jl)


The true profile shouldn't be incredibly spiky , and thus A\b should not induce many negative numbers in the A\b inversion). Thus the `active set` (numbers set to 0 because they violated the nonnegativity constraint) is usually small, so we usually don't need many iterations to find a solution

Because A is usually tridiagonal (maybe septadiagonal at worse for cubic splines) the inversions which might ordinarily be expensive are actually quite fast -- caching doesn't seem to improve performance.


Some sample timings with A = (1000x1000) and tridiagonal (but full matrix type, not LinearAlgebra.Tridiagonal) and M = (1000,) is fully dense and not sparse (few zeros):


@btime nonneg_lsq($A, $M; alg=:nnls) evals=2 samples=2  # doesn't support `tol`
    1.027 s (13 allocations: 7.68 MiB)

@btime nonneg_lsq($A, $M; alg=:nnls) evals=2 samples=2 


@btime nonneg_lsq($A, $M; alg=:fnnls, tol=1e-12) evals=2 samples=2 
    6.200 s (22914 allocations: 3.71 GiB)
@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-8) evals=2 samples=2
    6.258 s (22914 allocations: 3.71 GiB)


@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-8) evals=2 samples=2
    33.881 ms (48 allocations: 38.25 MiB)     : This is 
@btime nonneg_lsq($A, $M; alg=:nnls, tol=1e-12) evals=2 samples=2
    35.031 ms (48 allocations: 38.25 MiB)

So, :pivot is over 150x faster, with about 100x fewer allocations than :fnnls, but does allocate more than :mmls

For the pivot variants, :cache and :comb

@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-12, variant=:cache) evals=2 samples=2
    909.911 ms (67 allocations: 76.54 MiB)

@btime nonneg_lsq($A, $M; alg=:pivot, tol=1e-12, variant=:comb) evals=2 samples=2
    932.510 ms (113 allocations: 76.54 MiB)

I thk the active sets are just too small for caching things to be worth it


The only quirk is that :pivot can leave some very small negative numbers in the solution, which we have to resolve to 0, but that's very fast to do.

This probably also explains why :pivot was the default in the original code.

Note if you're doing fully random matrices for A, :pivot was slower, but then that's not really a realistic representation of what the workflow would be for regridding input variables like this where A is sparse and has a real structure from the grid.
=#

const default_conservative_interp_kwargs = (;
    preserve_monotonicity = true,
    enforce_positivity = false,
    nnls_alg = :pivot, # should be optimal for most cases, see the discussion above.
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
