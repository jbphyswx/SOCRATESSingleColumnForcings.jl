using Documenter: Documenter
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings, SOCRATESSingleColumnForcings as SSCF

#! format: off
api = Any[
    "Forcings" => "api/forcings.md",
    "Data I/O" => "api/io.md",
    "Interpolation" => "api/interpolation.md",
    "Thermodynamics" => "api/thermodynamics.md",
]

pages = Any[
    "Home" => "index.md",
    "Getting started" => "getting-started.md",
    "Column forcings" => "forcings.md",
    "Interpolation" => "interpolation.md",
    "Data and artifacts" => "data-and-artifacts.md",
    "API" => api,
    "References" => "references.md",
]

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    collapselevel = 2,
)
#! format: on

Documenter.makedocs(
    sitename = "SOCRATESSingleColumnForcings.jl",
    format = format,
    clean = true,
    modules = [SSCF, SSCF.Interpolation],
    pagesonly = true,
    # The package deliberately exports nothing (qualified access only), so many internal helpers
    # carry docstrings that are intentionally NOT surfaced on the public API pages. Keep the full
    # `:all` docstring-coverage check enabled (so genuine gaps still surface as warnings) but do not
    # fail the build over those internal-only docstrings.
    checkdocs = :all,
    warnonly = [:missing_docs],
    pages = pages,
)

Documenter.deploydocs(
    repo = "github.com/jbphyswx/SOCRATESSingleColumnForcings.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
)
