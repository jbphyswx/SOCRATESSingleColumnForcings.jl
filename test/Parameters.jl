"""
# TODO: replace with TurbulenceConvectionParameters from TurbulenceConvection
"""
module Parameters

# ADAPTED FROM  https://github.com/CliMA/TurbulenceConvection.jl/blob/99d35cafd158d8793fa744f773d5f942d6e74488/driver/parameter_set.jl
# and https://github.com/CliMA/TurbulenceConvection.jl/blob/7b4666baca418b00bb60e929de96fcc06100de57/src/Parameters.jl
# using SurfaceFluxes: SurfaceFluxes as SF
# using CloudMicrophysics: CloudMicrophysics as CM

abstract type AbstractTurbulenceConvectionParameters end
const ATCP = AbstractTurbulenceConvectionParameters


using CLIMAParameters: CLIMAParameters as CP
using Thermodynamics: Thermodynamics as TD

pairs = NamedTuple()

# toml_dict_default::CP.AbstractTOMLDict
FT = Float64
FTD = Float64

toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; param_pairs...)
TP = typeof(thermo_params)


aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
aliases = setdiff(aliases, ["thermo_params"])
pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
microphys_params = CM.Parameters.CloudMicrophysicsParameters{FTD, TP}(; pairs..., thermo_params)
MP = typeof(microphys_params)


# #####
# ##### TurbulenceConvection parameters
# #####

Base.@kwdef struct TurbulenceConvectionParameters{MP} <: ATCP
    # Omega::FT
    # planet_radius::FT
    # microph_scaling::FT
    # microph_scaling_dep_sub::FT
    # microph_scaling_melt::FT
    # microph_scaling_acnv::FT
    # microph_scaling_accr::FT
    microphys_params::MP
    # surf_flux_params::SFP
    # user_args::NamedTuple # not sure if this is completely necessary yet
end

thermodynamics_params(ps::Union{ATCP}) = CM.Parameters.thermodynamics_params(ps.microphys_params)
# surface_fluxes_params(ps::ATCP) = ps.surf_flux_params
# microphysics_params(ps::ATCP) = ps.microphys_params

# Base.eltype(::TurbulenceConvectionParameters{FT}) where {FT} = FT
# Omega(ps::ATCP) = ps.Omega
# planet_radius(ps::ATCP) = ps.planet_radius
# # TODO - microph_scaling is the factor for adjusting evaporation.
# # The name will be fixed in CLIMAParameters first.
# microph_scaling(ps::ATCP) = ps.microph_scaling
# microph_scaling_dep_sub(ps::ATCP) = ps.microph_scaling_dep_sub
# microph_scaling_melt(ps::ATCP) = ps.microph_scaling_melt
# microph_scaling_acnv(ps::ATCP) = ps.microph_scaling_acnv
# microph_scaling_accr(ps::ATCP) = ps.microph_scaling_accr

# #####
# ##### Forwarding parameters
# #####
# # Gonna try to figure out all the parameters we need to make a parameter set for what is in this module w/o too many external dependinces esp TC.jl, SF.jl, CM.jl, etc

for var in fieldnames(TD.Parameters.ThermodynamicsParameters)
    @eval $var(ps::ATCP) = TD.Parameters.$var(thermodynamics_params(ps))
end

# derived parameters
molmass_ratio(ps::Union{ATCP}) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
R_d(ps::Union{ATCP}) = TD.Parameters.R_d(thermodynamics_params(ps))
R_v(ps::Union{ATCP}) = TD.Parameters.R_v(thermodynamics_params(ps))
cp_d(ps::Union{ATCP}) = TD.Parameters.cp_d(thermodynamics_params(ps))
cv_v(ps::Union{ATCP}) = TD.Parameters.cv_v(thermodynamics_params(ps))
cv_l(ps::Union{ATCP}) = TD.Parameters.cv_l(thermodynamics_params(ps))

# ##### Forwarding SurfaceFluxes.jl

# # von_karman_const(ps::ATCP) = SF.Parameters.von_karman_const(surface_fluxes_params(ps))

# ##### Forwarding CloudMicrophysics.jl

ρ_cloud_liq(ps::ATCP) = CM.Parameters.ρ_cloud_liq(microphysics_params(ps))

param_set = TurbulenceConvectionParameters{MP}(; microphys_params)


end
