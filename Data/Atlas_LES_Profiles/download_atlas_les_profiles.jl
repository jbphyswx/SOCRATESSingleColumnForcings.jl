#=
Download LES forcing data from https://atmos.uw.edu/~ratlas/SOCRATES-LES-cases.html
=#

using NCDatasets: NCDatasets as NC  
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF


thisdir = @__DIR__ # doesn't seem to work to use @__DIR__ directly as a variable
const package_data_dir = dirname(thisdir)

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

# Artifact names: one per (flight, forcing) for the input & output forcing `.nc` files, plus a single
# shared metadata artifact holding `SOCRATES_summary.nc` + the two unique level grids. Each is declared
# as a lazy download artifact in the package `Artifacts.toml` (no data is duplicated across artifacts).
_forcing_tag(::ObsForcing) = "obs"
_forcing_tag(::ERA5Forcing) = "era5"
_inputs_artifact_name(flight::Int, ft::AbstractForcingType) =
    "atlas_les_inputs_" * lowercase(_rf_num(flight)) * "_" * _forcing_tag(ft) * "_v1"
_outputs_artifact_name(flight::Int, ft::AbstractForcingType) =
    "atlas_les_outputs_" * lowercase(_rf_num(flight)) * "_" * _forcing_tag(ft) * "_v1"
const _metadata_artifact_name = "atlas_les_metadata_v1"

# Resolve a declared (lazy) artifact by name to its on-disk path, downloading on first use from the
# `[[<name>.download]]` mirror(s) in `Artifacts.toml`. `LazyArtifacts` (loaded by the module) enables the
# on-demand download; nothing is written to the package's `Artifacts.toml` at runtime, so this works in
# read-only depots and offline once cached.
function _artifact_root(name::AbstractString)
    Pkg.Artifacts.ensure_artifact_installed(name, SSCF.artifacts_toml)
    return Pkg.Artifacts.artifact_path(Pkg.Artifacts.artifact_hash(name, SSCF.artifacts_toml))
end


# ---------------------------------------------------------------------------------------------------
# Runtime resolution: each artifact is a DECLARED lazy artifact in `Artifacts.toml` (git-tree-sha1 +
# `[[download]]` mirror(s) — GitHub release / Zenodo). `_artifact_root` downloads it on first use and
# returns its path. No `create_artifact`/`bind_artifact!` at runtime (that path — which mutated the
# package `Artifacts.toml` and broke in read-only depots / CI — is gone). Regenerating the tarballs from
# the raw Box files is a separate maintenance tool below (`regenerate_artifact_tarballs`).
# ---------------------------------------------------------------------------------------------------

"""
    atlas_les_inputs_root(flight, forcing_type)

On-disk directory of the input-forcing artifact for `(flight, forcing_type)` — one `.nc` at its root
(`RFNN_<obs|ERA5>…_SAM_input…nc`). Lazily downloaded on first use (see `Artifacts.toml`).
"""
atlas_les_inputs_root(flight::Int, forcing_type::AbstractForcingType) =
    _artifact_root(_inputs_artifact_name(flight, forcing_type))

"""
    atlas_les_outputs_root(flight, forcing_type)

On-disk directory of the output-forcing (LES) artifact for `(flight, forcing_type)` — one `.nc` at its
root. Lazily downloaded on first use.
"""
atlas_les_outputs_root(flight::Int, forcing_type::AbstractForcingType) =
    _artifact_root(_outputs_artifact_name(flight, forcing_type))

"""
    atlas_les_metadata_root()

On-disk directory of the shared metadata artifact: `SOCRATES_summary.nc` + the two level grids
(`192level-grd.txt`, `320level-grd.txt`). A flight's grid is `\$(grid_heights[flight])level-grd.txt`.
"""
atlas_les_metadata_root() = _artifact_root(_metadata_artifact_name)

"Path to `SOCRATES_summary.nc` (in the shared metadata artifact). `flight_number` is accepted for call-site compatibility; the summary is flight-independent."
atlas_socrates_summary_file(flight_number::Int) = joinpath(atlas_les_metadata_root(), "SOCRATES_summary.nc")

"Path to the level-grid file for `flight` (in the shared metadata artifact), selected via `grid_heights`."
atlas_grid_file(flight::Int) = joinpath(atlas_les_metadata_root(), string(SSCF.grid_heights[flight]) * "level-grd.txt")

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




# ---------------------------------------------------------------------------------------------------
# Raw retrieval / regeneration from Box (NOT the runtime path). The package serves data at runtime via
# the declared lazy artifacts in `Artifacts.toml`; these helpers fetch the raw per-file data from the
# canonical Box source into a directory so the artifact tarballs can be (re)built and re-uploaded:
#   dir = download_atlas_les_inputs(mktempdir()); download_atlas_les_outputs(dir)
#   # then per (flight,forcing): Pkg.Artifacts.create_artifact() do d; cp(nc, ...); end; archive_artifact(h, "name.tar.gz")
# The Box links + `_download_first` are the retained canonical raw source.
# ---------------------------------------------------------------------------------------------------

"""
    download_atlas_les_inputs(destdir; flight_numbers, forcing_types = (:obs_data, :ERA5_data))

Download the raw Atlas LES *input* files (per forcing `.nc` + level grids) from Box into `destdir`.
"""
function download_atlas_les_inputs(
    destdir::AbstractString;
    flight_numbers::Union{AbstractArray{Int}, Tuple{Vararg{Int}}} = flight_numbers,
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data, :ERA5_data),
    SOCRATES_LES_inputs_Box_links::Dict{String, String} = SOCRATES_LES_inputs_Box_links,
)
    mkpath(destdir)
    for flight in flight_numbers
        RF_num = _rf_num(flight)
        for forcing_type in forcing_types
            fn = forcing_type === :obs_data ? RF_num * "_obs-based_SAM_input.nc" :
                 forcing_type === :ERA5_data ? RF_num * "_ERA5-based_SAM_input_mar18_2022.nc" :
                 error("forcing_type must be :obs_data or :ERA5_data")
            urls = (get(SOCRATES_LES_inputs_Box_links, fn, nothing), uw_atlas_base_url * fn)
            _download_first(urls, joinpath(destdir, fn)) || @warn "input not found on Box/UW: $fn"
        end
        gfn = string(SSCF.grid_heights[flight]) * "level-grd.txt"
        urls = (get(SOCRATES_LES_inputs_Box_links, gfn, nothing), uw_atlas_base_url * gfn)
        _download_first(urls, joinpath(destdir, gfn)) || @warn "grid not found on Box/UW: $gfn"
    end
    return destdir
end

"""
    download_atlas_les_outputs(destdir; flight_numbers, forcing_types = (:obs_data, :ERA5_data))

Download the raw Atlas LES *output* files from Box into `destdir`.
"""
function download_atlas_les_outputs(
    destdir::AbstractString;
    flight_numbers::Union{AbstractArray{Int}, Tuple{Vararg{Int}}} = flight_numbers,
    forcing_types::Union{AbstractArray{Symbol}, Tuple{Vararg{Symbol}}} = (:obs_data, :ERA5_data),
    SOCRATES_LES_outputs_Box_links::Dict{String, String} = SOCRATES_LES_outputs_Box_links,
)
    mkpath(destdir)
    suffix = "_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc"
    for flight in flight_numbers
        RF_num = _rf_num(flight)
        for forcing_type in forcing_types
            tag = forcing_type === :obs_data ? "obs" : forcing_type === :ERA5_data ? "ERA5" :
                  error("forcing_type must be :obs_data or :ERA5_data")
            fn = RF_num * "_" * tag * suffix
            urls = (get(SOCRATES_LES_outputs_Box_links, fn, nothing), uw_atlas_base_url * RF_num * "_output/" * lowercase(tag) * "/SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc")
            _download_first(urls, joinpath(destdir, fn)) || @warn "output not found on Box/UW: $fn"
        end
    end
    return destdir
end
