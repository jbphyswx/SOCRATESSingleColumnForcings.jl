# Documentation (Documenter.jl)

User guides are built with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).

## Build locally

```bash
cd SOCRATESSingleColumnForcings.jl
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` in a browser.

## Source layout

```
docs/
  make.jl
  src/
    index.md
    getting-started.md
    forcings.md
    interpolation.md
    data-and-artifacts.md
    references.md
```

## Deploy

CI (`.github/workflows/docs.yml`) builds and deploys to GitHub Pages on pushes to `main` and on tags. Requires a `DOCUMENTER_KEY` repository secret for SSH deploy.

Deployed site: [jbphyswx.github.io/SOCRATESSingleColumnForcings.jl](https://jbphyswx.github.io/SOCRATESSingleColumnForcings.jl/dev/)
