#=
Download LES forcing data from https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html
=#

import NCDatasets as NC
import SOCRATESSingleColumnForcings as SSCF


thisdir = @__DIR__ # doesn't seem to work to use @__DIR__ directly as a variable
const package_data_dir = dirname(thisdir)

function _copytree_contents!(src::AbstractString, dst::AbstractString)
    for name in readdir(src)
        src_path = joinpath(src, name)
        dst_path = joinpath(dst, name)
        if isdir(src_path)
            mkpath(dst_path)
            _copytree_contents!(src_path, dst_path)
        else
            cp(src_path, dst_path; force = true)
        end
    end
    return dst
end

function _download_first(urls, savepath::AbstractString)
    found = false
    for url in urls
        if isnothing(url)
            continue
        end
        try
            Downloads.download(url, savepath)
            found = true
            @info "Found $(url) for $savepath"
            break
        catch
            @warn "Did not find $(url)"
        end
    end
    return found
end

_rf_num(flight::Int) = "RF" * string(flight, pad = 2)
_inputs_artifact_name(flight::Int) = "atlas_les_inputs_" * lowercase(_rf_num(flight)) * "_v1"
_outputs_artifact_name(flight::Int) = "atlas_les_outputs_" * lowercase(_rf_num(flight)) * "_v1"

function _artifact_dir(name::String)
    hash = Artifacts.artifact_hash(name, SSCF.artifacts_toml)
    return isnothing(hash) ? nothing : Artifacts.artifact_path(hash)
end

function _required_input_files(flight::Int, forcing_types)
    RF_num = _rf_num(flight)
    req = [RF_num * "_grd.txt"]
    for forcing_type in forcing_types
        if forcing_type == :obs_data
            push!(req, RF_num * "_obs-based_SAM_input.nc")
        elseif forcing_type == :ERA5_data
            push!(req, RF_num * "_ERA5-based_SAM_input_mar18_2022.nc")
        else
            error("forcing_type must be :obs_data or :ERA5_data")
        end
    end
    return req
end

function _required_output_files(flight::Int, forcing_types)
    RF_num = _rf_num(flight)
    req = [RF_num * "_grd.txt"]
    for forcing_type in forcing_types
        if forcing_type == :obs_data
            push!(req, RF_num * "_Obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc")
        elseif forcing_type == :ERA5_data
            push!(req, RF_num * "_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc")
        else
            error("forcing_type must be :obs_data or :ERA5_data")
        end
    end
    return req
end

function atlas_les_inputs_root(
    flight::Int;
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data,),
    SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links,
)
    name = _inputs_artifact_name(flight)
    existing = _artifact_dir(name)
    input_data_dir = isnothing(existing) ? nothing : joinpath(existing, "Input_Data")
    required = _required_input_files(flight, forcing_types)
    summary_artifact_file = isnothing(existing) ? nothing : joinpath(existing, "SOCRATES_summary.nc")
    has_required = !isnothing(existing) &&
                   all(isfile(joinpath(input_data_dir, file)) for file in required) &&
                   isfile(summary_artifact_file)
    if has_required
        return existing
    end

    RF_num = _rf_num(flight)
    new_hash = Pkg.Artifacts.create_artifact() do artifact_dir
        if !isnothing(existing)
            _copytree_contents!(existing, artifact_dir)
        end

        out_dir = joinpath(artifact_dir, "Input_Data")
        mkpath(out_dir)

        for forcing_type in forcing_types
            if forcing_type == :obs_data
                obs_filename = RF_num * "_obs-based_SAM_input.nc"
                obs_savepath = joinpath(out_dir, obs_filename)
                if !isfile(obs_savepath)
                    urls = (
                        get(SOCRATES_LES_inputs_Box_links, obs_filename, nothing),
                        uw_atlas_base_url * obs_filename,
                    )
                    found = _download_first(urls, obs_savepath)
                    !found && @error "Did not find obs forcing for $RF_num in $(urls)"
                end
            elseif forcing_type == :ERA5_data
                ERA5_filename = RF_num * "_ERA5-based_SAM_input_mar18_2022.nc"
                ERA5_savepath = joinpath(out_dir, ERA5_filename)
                if !isfile(ERA5_savepath)
                    urls = (
                        get(SOCRATES_LES_inputs_Box_links, ERA5_filename, nothing),
                        uw_atlas_base_url * ERA5_filename,
                    )
                    found = _download_first(urls, ERA5_savepath)
                    !found && @error "Did not find ERA5 forcing for $RF_num in $(urls)"
                end
            else
                error("forcing_type must be :obs_data or :ERA5_data")
            end
        end

        grid_filename = RF_num * "_grd.txt"
        grid_savepath = joinpath(out_dir, grid_filename)
        if !isfile(grid_savepath)
            grid_height = SSCF.grid_heights[flight]
            urls = (
                get(SOCRATES_LES_inputs_Box_links, string(grid_height) * "level-grd.txt", nothing),
                uw_atlas_base_url * string(grid_height) * "level-grd.txt",
            )
            found = _download_first(urls, grid_savepath)
            !found && @error "Did not find grid file for $RF_num in $(urls)"
        end

        summary_src = joinpath(package_data_dir, "SOCRATES_summary.nc")
        summary_dst = joinpath(artifact_dir, "SOCRATES_summary.nc")
        if isfile(summary_src) && !isfile(summary_dst)
            cp(summary_src, summary_dst; force = true)
        end
    end

    Pkg.Artifacts.bind_artifact!(SSCF.artifacts_toml, name, new_hash; force = true)
    return Artifacts.artifact_path(new_hash)
end

function atlas_les_outputs_root(
    flight::Int;
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data, :ERA5_data),
    SOCRATES_LES_outputs_Box_links::Dict{String, String} = SOCRATES_LES_outputs_Box_links,
    SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links,
)
    name = _outputs_artifact_name(flight)
    existing = _artifact_dir(name)
    output_data_dir = isnothing(existing) ? nothing : joinpath(existing, "Output_Data")
    required = _required_output_files(flight, forcing_types)
    has_required = !isnothing(existing) && all(isfile(joinpath(output_data_dir, file)) for file in required)
    if has_required
        return existing
    end

    RF_num = _rf_num(flight)
    new_hash = Pkg.Artifacts.create_artifact() do artifact_dir
        if !isnothing(existing)
            _copytree_contents!(existing, artifact_dir)
        end

        out_dir = joinpath(artifact_dir, "Output_Data")
        mkpath(out_dir)

        for forcing_type in forcing_types
            if forcing_type == :obs_data
                obs_filename = RF_num * "_output/obs/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"
                obs_savepath = joinpath(out_dir, RF_num * "_Obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc")
                if !isfile(obs_savepath)
                    urls = (
                        get(
                            SOCRATES_LES_outputs_Box_links,
                            RF_num * "_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc",
                            nothing,
                        ),
                        uw_atlas_base_url * obs_filename,
                    )
                    found = _download_first(urls, obs_savepath)
                    !found && @error "Did not find obs forcing for $RF_num in $(urls)"
                end
            elseif forcing_type == :ERA5_data
                ERA5_filename = RF_num * "_output/era5/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"
                ERA5_savepath = joinpath(out_dir, RF_num * "_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc")
                if !isfile(ERA5_savepath)
                    urls = (
                        get(
                            SOCRATES_LES_outputs_Box_links,
                            RF_num * "_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc",
                            nothing,
                        ),
                        uw_atlas_base_url * ERA5_filename,
                    )
                    found = _download_first(urls, ERA5_savepath)
                    !found && @error "Did not find ERA5 forcing for $RF_num in either $(urls)"
                end
            else
                error("forcing_type must be :obs_data or :ERA5_data")
            end
        end

        grid_filename = RF_num * "_grd.txt"
        grid_savepath = joinpath(out_dir, grid_filename)
        if !isfile(grid_savepath)
            grid_height = SSCF.grid_heights[flight]
            urls = (
                get(SOCRATES_LES_inputs_Box_links, string(grid_height) * "level-grd.txt", nothing),
                uw_atlas_base_url * string(grid_height) * "level-grd.txt",
            )
            found = _download_first(urls, grid_savepath)
            !found && @error "Did not find output-grid file for $RF_num in $(urls)"
        end
    end

    Pkg.Artifacts.bind_artifact!(SSCF.artifacts_toml, name, new_hash; force = true)
    return Artifacts.artifact_path(new_hash)
end

function atlas_socrates_summary_file(flight_number::Int)
    artifact_dir = atlas_les_inputs_root(flight_number; forcing_types = (:obs_data,))
    summary_file = joinpath(artifact_dir, "SOCRATES_summary.nc")
    if isfile(summary_file)
        return summary_file
    end
    return joinpath(package_data_dir, "SOCRATES_summary.nc")
end

# download forcings for each flight

const uw_atlas_base_url = "https://atmos.uw.edu/~ratlas/"


const SOCRATES_flight_observations_Box_links = Dict{Int, String}( # The raw observational data from the SOCRATES flights, for ready access on any machine...
    # They're big files so not storing them in the git repo...
    1 => "https://caltech.box.com/shared/static/4sh03gadjq6acary69qamgjbmynfio2f.nc",
    2 => "https://caltech.box.com/shared/static/urhtegy7dccnp8hav4my90kzrr3e4wmt.nc",
    3 => "https://caltech.box.com/shared/static/xp4b4p2ef523bqmzuejbecmruu8daj6s.nc",
    4 => "https://caltech.box.com/shared/static/f41qh4jjlnvs08y676hhq60epcta2oip.nc",
    5 => "https://caltech.box.com/shared/static/ty0h3mxrj6myyun3qstfoj1ef0o8009m.nc",
    6 => "https://caltech.box.com/shared/static/e3yk43xfwyxyix0yqsh53utpmjynzkzr.nc",
    7 => "https://caltech.box.com/shared/static/34bz9z9hlolqwmer02g2qz4uu13qn13s.nc",
    8 => "https://caltech.box.com/shared/static/pua0gf97772xn5cdq256bztriu4q5kcz.nc",
    9 => "https://caltech.box.com/shared/static/1vdwhg0oiwtzd21vqwzt5sm6zsh8yng9.nc",
    10 => "https://caltech.box.com/shared/static/cw40hwxoegcmjxksja4gxqfpa22mwj9p.nc",
    11 => "https://caltech.box.com/shared/static/cz7qhv4ub5kbsxloxydf2kyykd9e5xdn.nc",
    12 => "https://caltech.box.com/shared/static/mugnp8iqfcccksxd8uojqbufcn8waexh.nc",
    13 => "https://caltech.box.com/shared/static/qzj14slfrknboi4iflen8t5gqfgmi595.nc",
    14 => "https://caltech.box.com/shared/static/kxcaogfc1m5eoopb5z3gedffv8gjizax.nc",
    15 => "https://caltech.box.com/shared/static/pfp8egkqcxhbyi1bxqfgk28g7neih7d8.nc",
)


# We have alias links for each and each file has a direct link that's accessible to programs w/o redirect
const SOCRATES_LES_inputs_Box_links = Dict{String, String}( # We have to save each link individually bc there's no way to traverse the box subfolders...
    # "192level-grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-192-level-grd",
    # "320level-grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-320-level-grd",
    "192level-grd.txt" => "https://caltech.box.com/shared/static/z3bys5xvcycgufe5z5iangf0rynn19fz.txt",
    "320level-grd.txt" => "https://caltech.box.com/shared/static/jz97xm12z3k1x1g7ypmo49me4c6v5abd.txt",
    #
    # "RF01_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF01-grd",
    # "RF09_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF09-grd",
    # "RF10_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF10-grd",
    # "RF11_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF11-grd",
    # "RF12_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF12-grd",
    # "RF13_grd.txt" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF13-grd",
    "RF01_grd.txt" => "https://caltech.box.com/shared/static/36dma4g1zvrgr8p9kzdneax6c256zm12.txt",
    "RF09_grd.txt" => "https://caltech.box.com/shared/static/cqztxdezusdi60joz8d3cnb47f9d63xq.txt",
    "RF10_grd.txt" => "https://caltech.box.com/shared/static/f9flfzv02opybzg7bsoi0vmwgvw8hh19.txt",
    "RF11_grd.txt" => "https://caltech.box.com/shared/static/5opm2r7tbix44zp74ohs9g5b620x5c4q.txt",
    "RF12_grd.txt" => "https://caltech.box.com/shared/static/ds4yq4ehhjchhiknsvwfbq0n1eopnclu.txt",
    "RF13_grd.txt" => "https://caltech.box.com/shared/static/p6yxgmu36ys9s95f9xr0xgbghyj9dvmd.txt",

    #
    # "RF01_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-jbphyswx-SSCFjl-RF01-ERA-SAMin",
    # "RF01_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF01-obs-SAMin",
    "RF01_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/yhfsgrg5a5g5yufpomk9oj4qqoha1ky8.nc",
    "RF01_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/qy40xltlvf2j8co5b5i9iwvgvd5nwmbq.nc",
    #
    # "RF09_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF09-ERA-SAMin",
    # "RF09_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF09-obs-SAMin",
    "RF09_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/4lpiceoom3i2pi1ea7b1znql2701q8gq.nc",
    "RF09_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/vr3src4yq06p1fwbd2jrm5uxfpj9r8z1.nc",
    #
    # "RF10_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF10-ERA-SAMin",
    # "RF10_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF10-obs-SAMin",
    "RF10_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/fwtkcsufgeesqce0hfa6ddamx93o5ti5.nc",
    "RF10_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/gce9r1q315gvze4216h6nx2ntlz33n07.nc",
    #
    # "RF11_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF11-ERA-SAMin",
    # "RF11_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF11-obs-SAMMin",
    "RF11_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/3hmn0nvzj9immd3ssa4gnbbreazg56jx.nc",
    # "RF11_obs-based_SAM_input.nc" => "",
    #
    # "RF12_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF12-ERA-SAMin",
    # "RF12_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF12-obs-SAMMin",
    "RF12_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/zjhihcxk5g7uq8hrnzeyxckurncjtifb.nc",
    "RF12_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/h9dcw7at13k154v800s57ubzuivkbjwy.nc",
    #
    # "RF13_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF13-ERA-SAMin",
    # "RF13_obs-based_SAM_input.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF13-obs-SAMMin",
    "RF13_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/hdy2x2cis97bs79f8p2qzniwddl6a84d.nc",
    "RF13_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/cczvc3zjb93vhvy7b5tm01y1l01s7hni.nc",
)

const SOCRATES_LES_outputs_Box_links = Dict{String, String}( # We have to save each link individually bc there's no way to traverse the box subfolders...
    # "RF01_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF01-ERA-SAMou",
    # "RF01_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF01-obs-SAMou",
    "RF01_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/x3qsjkglgshfb5792r8ra7idk9g9g9ju.nc",
    "RF01_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/lm6iqy26ey9xce5wke6kqggc5gd6tv2c.nc",
    #
    # "RF09_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF09-ERA-SAMou",
    # "RF09_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF09-obs-SAMou",
    "RF09_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/ezw97njysbj7j42tij44zzn3ltqylwes.nc",
    "RF09_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/eothh18t47lgu8cru8asab2295mlk524.nc",
    #
    # "RF10_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF10-ERA-SAMou",
    # "RF10_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF10-obs-SAMou",
    "RF10_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/xv8arlccio3jpzks5bpvx2mkxeqwb6t6.nc",
    "RF10_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/0niw6brvrcxr1mk8go2bkzvfw0nw3ifw.nc",
    #
    # "RF11_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF11-ERA-SAMou",
    # "RF11_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF11-obs-SAMou",
    "RF11_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/d152r6uo78kkblxq9z2hz01skgs41tbi.nc",
    # "RF11_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "",
    #
    # "RF12_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF12-ERA-SAMou",
    # "RF12_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF12-obs-SAMou",
    "RF12_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/cijfzktn35ja2ifeboosxbbn1p92gmoo.nc",
    "RF12_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/rxtphfcu9ay0ceu2htntognpjdzpcgdr.nc",
    #
    # "RF13_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF13-ERA-SAMou",
    # "RF13_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/v/jbphyswx-SSCFjl-RF13-obs-SAMou",
    "RF13_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/eyftsdsc0ux5ulx4nsnjp38v1h3wjqab.nc",
    "RF13_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/z7gwogfk5uvmrqotd8br2yqjib1ydlbc.nc",
)




"""
    download_atlas_les_inputs(; flight_numbers::AbstractArray{Int} = flight_numbers, forcing_type::Symbol = :obs_data, SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links)

"""
function download_atlas_les_inputs(;
    flight_numbers::Union{AbstractArray{Int}, Tuple{Vararg{Int}}} = flight_numbers,
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data,),
    SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links,
)
    for flight in flight_numbers
        atlas_les_inputs_root(
            flight;
            forcing_types = forcing_types,
            SOCRATES_LES_inputs_Box_links = SOCRATES_LES_inputs_Box_links,
        )
    end
    return nothing
end






"""
    download_atlas_les_profiles(; flight_numbers::AbstractArray{Int} = flight_numbers, forcing_types::AbstractArray{Symbol} = [:obs_data, :ERA5_data], SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links, SOCRATES_LES_outputs_Box_links::Dict{String, String} = SOCRATES_LES_outputs_Box_links)

    These all have the same fiilename so we have to manually append the flight number to the filename

"""

function download_atlas_les_outputs(;
    flight_numbers::Union{AbstractArray{Int}, Tuple{Vararg{Int}}} = flight_numbers,
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data, :ERA5_data),
    SOCRATES_LES_outputs_Box_links::Dict{String, String} = SOCRATES_LES_outputs_Box_links,
)
    for flight in flight_numbers
        atlas_les_outputs_root(
            flight;
            forcing_types = forcing_types,
            SOCRATES_LES_outputs_Box_links = SOCRATES_LES_outputs_Box_links,
            SOCRATES_LES_inputs_Box_links = SOCRATES_LES_inputs_Box_links,
        )
    end
    return nothing
end
