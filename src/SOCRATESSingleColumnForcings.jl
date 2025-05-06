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

# include our files
const FT = Float64
const flight_numbers = [1, 9, 10, 11, 12, 13]
const forcing_types = [:obs_data, :ERA5_data] # maybe change these to [:obs,:ERA5] later? would need to mirror in Cases.jl in TC.jl
const TDPS = TD.Parameters.ThermodynamicsParameters
const TDTS = TD.ThermodynamicState

# Two cases with shallow cloud-topped boundary layers, RF12 and RF13, are run on a 192-level vertical grid.
# The other four cases have clouds extending through deeper boundary layers; they are run on a 320-level vertical grid.
const grid_heights = Dict(1 => 320, 9 => 320, 10 => 320, 11 => 320, 12 => 192, 13 => 192)

include(joinpath("..", "Data", "Atlas_LES_Profiles", "download_atlas_les_profiles.jl")) # not in src so don't include automatically
include("open_atlas_les_inputs.jl")
include("open_atlas_les_outputs.jl")
include("helper_functions.jl")
include("process_case.jl")
# include("interpolating_methods.jl") # currently only helper_functions.jl reads these
# include("les_reader_helper.jl") # currently only helper_functions.jl reads these
include("get_LES_reference_profiles.jl")

end
