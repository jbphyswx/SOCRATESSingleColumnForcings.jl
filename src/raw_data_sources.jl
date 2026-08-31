#=
Raw upstream retrieval from the canonical Box / UW sources. This is NOT the runtime path — the
package serves data through the lazy artifacts resolved in `artifacts.jl`. These helpers fetch the
raw per-file data into a directory so the artifact tarballs can be rebuilt and re-uploaded:

    dir = download_atlas_les_inputs(mktempdir()); download_atlas_les_outputs(dir)
    # then per (flight, forcing): Pkg.Artifacts.create_artifact() do d; cp(nc, ...); end
    #                             archive_artifact(h, "name.tar.gz")
=#

const uw_atlas_base_url = "https://atmos.uw.edu/~ratlas/"

function _download_first(urls, savepath::AbstractString)
    for url in urls
        isnothing(url) && continue
        try
            Downloads.download(url, savepath)
            @info "Found $(url) for $savepath"
            return true
        catch
            @warn "Did not find $(url)"
        end
    end
    return false
end

# We have alias links for each and each file has a direct link that's accessible to programs w/o redirect
const SOCRATES_LES_inputs_Box_links = Dict{String, String}( # We have to save each link individually bc there's no way to traverse the box subfolders...
    "192level-grd.txt" => "https://caltech.box.com/shared/static/z3bys5xvcycgufe5z5iangf0rynn19fz.txt",
    "320level-grd.txt" => "https://caltech.box.com/shared/static/jz97xm12z3k1x1g7ypmo49me4c6v5abd.txt",
    #
    "RF01_grd.txt" => "https://caltech.box.com/shared/static/36dma4g1zvrgr8p9kzdneax6c256zm12.txt",
    "RF09_grd.txt" => "https://caltech.box.com/shared/static/cqztxdezusdi60joz8d3cnb47f9d63xq.txt",
    "RF10_grd.txt" => "https://caltech.box.com/shared/static/f9flfzv02opybzg7bsoi0vmwgvw8hh19.txt",
    "RF11_grd.txt" => "https://caltech.box.com/shared/static/5opm2r7tbix44zp74ohs9g5b620x5c4q.txt",
    "RF12_grd.txt" => "https://caltech.box.com/shared/static/ds4yq4ehhjchhiknsvwfbq0n1eopnclu.txt",
    "RF13_grd.txt" => "https://caltech.box.com/shared/static/p6yxgmu36ys9s95f9xr0xgbghyj9dvmd.txt",
    #
    "RF01_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/yhfsgrg5a5g5yufpomk9oj4qqoha1ky8.nc",
    "RF01_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/qy40xltlvf2j8co5b5i9iwvgvd5nwmbq.nc",
    #
    "RF09_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/4lpiceoom3i2pi1ea7b1znql2701q8gq.nc",
    "RF09_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/vr3src4yq06p1fwbd2jrm5uxfpj9r8z1.nc",
    #
    "RF10_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/fwtkcsufgeesqce0hfa6ddamx93o5ti5.nc",
    "RF10_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/gce9r1q315gvze4216h6nx2ntlz33n07.nc",
    #
    # RF11 is ERA5-only; there is no obs-based input file upstream.
    "RF11_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/3hmn0nvzj9immd3ssa4gnbbreazg56jx.nc",
    #
    "RF12_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/zjhihcxk5g7uq8hrnzeyxckurncjtifb.nc",
    "RF12_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/h9dcw7at13k154v800s57ubzuivkbjwy.nc",
    #
    "RF13_ERA5-based_SAM_input_mar18_2022.nc" => "https://caltech.box.com/shared/static/hdy2x2cis97bs79f8p2qzniwddl6a84d.nc",
    "RF13_obs-based_SAM_input.nc" => "https://caltech.box.com/shared/static/cczvc3zjb93vhvy7b5tm01y1l01s7hni.nc",
)

const SOCRATES_LES_outputs_Box_links = Dict{String, String}( # We have to save each link individually bc there's no way to traverse the box subfolders...
    "RF01_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/x3qsjkglgshfb5792r8ra7idk9g9g9ju.nc",
    "RF01_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/lm6iqy26ey9xce5wke6kqggc5gd6tv2c.nc",
    #
    "RF09_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/ezw97njysbj7j42tij44zzn3ltqylwes.nc",
    "RF09_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/eothh18t47lgu8cru8asab2295mlk524.nc",
    #
    "RF10_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/xv8arlccio3jpzks5bpvx2mkxeqwb6t6.nc",
    "RF10_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/0niw6brvrcxr1mk8go2bkzvfw0nw3ifw.nc",
    #
    # RF11 is ERA5-only; there is no obs-based output file upstream.
    "RF11_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/d152r6uo78kkblxq9z2hz01skgs41tbi.nc",
    #
    "RF12_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/cijfzktn35ja2ifeboosxbbn1p92gmoo.nc",
    "RF12_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/rxtphfcu9ay0ceu2htntognpjdzpcgdr.nc",
    #
    "RF13_ERA5_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/eyftsdsc0ux5ulx4nsnjp38v1h3wjqab.nc",
    "RF13_obs_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc" => "https://caltech.box.com/shared/static/z7gwogfk5uvmrqotd8br2yqjib1ydlbc.nc",
)

"""
    download_atlas_les_inputs(destdir; flight_numbers, forcing_types = forcing_types)

Download the raw Atlas LES *input* files (per forcing `.nc` + the selected flight grid) into `destdir`.
Maintenance tool for rebuilding the artifacts; the runtime path resolves through `Artifacts.toml`.
"""
function download_atlas_les_inputs(
    destdir::AbstractString;
    flight_numbers::Union{AbstractArray{<:Integer}, Tuple{Vararg{Integer}}} = flight_numbers,
    forcing_types::Union{AbstractArray{<:AbstractForcingType}, Tuple{Vararg{AbstractForcingType}}} = forcing_types,
    links::Dict{String, String} = SOCRATES_LES_inputs_Box_links,
)
    mkpath(destdir)
    for flight in flight_numbers
        RF_num = _rf_num(flight)
        for forcing_type in forcing_types
            fn = RF_num * "_" * _input_file_tag(forcing_type) * ".nc"
            urls = (get(links, fn, nothing), uw_atlas_base_url * fn)
            _download_first(urls, joinpath(destdir, fn)) || @warn "input not found on Box/UW: $fn"
        end
        gfn = _atlas_grid_filename(flight)
        urls = (get(links, gfn, nothing), uw_atlas_base_url * gfn)
        _download_first(urls, joinpath(destdir, gfn)) || @warn "grid not found on Box/UW: $gfn"
    end
    return destdir
end

"""
    download_atlas_les_outputs(destdir; flight_numbers, forcing_types = forcing_types)

Download the raw Atlas LES *output* files into `destdir`. Maintenance tool for rebuilding the
artifacts; the runtime path resolves through `Artifacts.toml`.
"""
function download_atlas_les_outputs(
    destdir::AbstractString;
    flight_numbers::Union{AbstractArray{<:Integer}, Tuple{Vararg{Integer}}} = flight_numbers,
    forcing_types::Union{AbstractArray{<:AbstractForcingType}, Tuple{Vararg{AbstractForcingType}}} = forcing_types,
    links::Dict{String, String} = SOCRATES_LES_outputs_Box_links,
)
    mkpath(destdir)
    for flight in flight_numbers
        RF_num = _rf_num(flight)
        for forcing_type in forcing_types
            tag = _output_file_tag(forcing_type)
            fn = RF_num * "_" * tag * _OUTPUT_FILE_SUFFIX
            urls = (
                get(links, fn, nothing),
                uw_atlas_base_url * RF_num * "_output/" * lowercase(tag) * "/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc",
            )
            _download_first(urls, joinpath(destdir, fn)) || @warn "output not found on Box/UW: $fn"
        end
    end
    return destdir
end
