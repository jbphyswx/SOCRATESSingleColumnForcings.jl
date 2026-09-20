#=
Read-boundary unit handling.

Every variable of the Atlas archive is declared here with the units the file carries, the units it is
converted to, and the transform between them.

This is needed because the files' `units` attributes are variously correct (`lev:units = "Pa"`),
absent (`Ps` carries an empty string), wrong (`divq` is labelled `K/s` while its `long_name` is a
moisture advection; `TQ` is labelled `K2` while it is a covariance of two different quantities),
self-referential (`SHFOBS:units = "SHFOBS"`), spelled several ways for one unit (`K/day` / `K/d`,
`mm/day` / `mm/d`), and — for the 82 microphysical process rates — a literal `?/kg/day` shared by
mass rates, number rates and inverse lengths alike. Each spec therefore records BOTH the raw units
expected on the file and the authority that fixed the classification when the attribute could not.
=#

# --- conversion primitives ------------------------------------------------------------------- #
# Scalar factors are written as functions

_div_1000(x::FT) where {FT} = x / FT(1000)

const g_to_kg = _div_1000
hPa_to_Pa(x::FT) where {FT} = x * FT(100)
perday_to_persec(x::FT) where {FT} = x / FT(86400)
percm3_to_perm3(x::FT) where {FT} = x * FT(1e6)
day_to_s(x::FT) where {FT} = x * FT(86400)
const mm_to_m = _div_1000
km_to_m(x::FT) where {FT} = x * FT(1000)
km2_to_m2(x::FT) where {FT} = x * FT(1e6)
micron_to_m(x::FT) where {FT} = x * FT(1e-6)
percent_to_fraction(x::FT) where {FT} = x / FT(100)
g2_to_kg2(x::FT) where {FT} = x / FT(1e6)

for f in (
    :g_to_kg, :hPa_to_Pa, :perday_to_persec, :percm3_to_perm3, :day_to_s,
    :km_to_m, :km2_to_m2, :micron_to_m, :percent_to_fraction, :g2_to_kg2,
)
    @eval $f(x::AbstractArray{FT}) where {FT <: AbstractFloat} = $f.(x)
end

"""
    wt_to_qt(wt)

Total-water mass mixing ratio [kg/kg dry air] → total-water specific humidity [kg/kg moist air],
`q_t = w_t / (1 + w_t)`.
"""
wt_to_qt(wt::FT) where {FT <: AbstractFloat} = wt / (one(FT) + wt)
wt_to_qt(wt::AbstractArray{FT}) where {FT <: AbstractFloat} = wt_to_qt.(wt)

"""
    wt_to_qt!(wt)

In-place [`wt_to_qt`](@ref).
"""
function wt_to_qt!(wt::AbstractArray)
    @inbounds for i in eachindex(wt)
        wt[i] = wt_to_qt(wt[i])
    end
    return wt
end

"""
    wc_to_qc(wc, qt)

Mass mixing ratio of a water species [kg/kg dry air] → its specific content [kg/kg moist air],
`q_x = w_x (1 - q_t)`, where `qt` is the TOTAL water specific humidity of the same air.

Holds for every species that is part of the total water, vapour included, since
`q_x = w_x/(1 + w_t)` and `1 - q_t = 1/(1 + w_t)`.
"""
wc_to_qc(wc::FT, qt::FT) where {FT <: AbstractFloat} = wc * (one(FT) - qt)
wc_to_qc(wc::AbstractArray{FT}, qt::AbstractArray{FT}) where {FT <: AbstractFloat} = wc_to_qc.(wc, qt)

"""
    qc_to_wc(qc, qt)

Inverse of [`wc_to_qc`](@ref), `w_x = q_x / (1 - q_t)`.
"""
qc_to_wc(qc::FT, qt::FT) where {FT <: AbstractFloat} = qc / (one(FT) - qt)
qc_to_wc(qc::AbstractArray{FT}, qt::AbstractArray{FT}) where {FT <: AbstractFloat} = qc_to_wc.(qc, qt)

"""
    wc2_to_qc2(wc2, qt)

Second moment of a species mixing ratio → second moment of its specific content: the conversion is
linear in `(1 - q_t)`, so its square takes `(1 - q_t)^2`.
"""
wc2_to_qc2(wc2::FT, qt::FT) where {FT <: AbstractFloat} = wc2 * (one(FT) - qt)^2
wc2_to_qc2(wc2::AbstractArray{FT}, qt::AbstractArray{FT}) where {FT <: AbstractFloat} = wc2_to_qc2.(wc2, qt)

"""
    wm2_to_kgm2s(F, L)

Water flux reported as an energy flux [W/m2] → mass flux [kg/m2/s], dividing by the latent heat `L`
[J/kg] of the phase change the reporting convention used.
"""
wm2_to_kgm2s(F::FT, L::FT) where {FT <: AbstractFloat} = F / L
wm2_to_kgm2s(F::AbstractArray{FT}, L::FT) where {FT <: AbstractFloat} = wm2_to_kgm2s.(F, L)

"""
    K_to_Jkg(x, cp)

Static energy stored divided by the dry-air heat capacity [K] → specific energy [J/kg].
"""
K_to_Jkg(x::FT, cp::FT) where {FT <: AbstractFloat} = x * cp
K_to_Jkg(x::AbstractArray{FT}, cp::FT) where {FT <: AbstractFloat} = K_to_Jkg.(x, cp)

"""
    permicron_to_perm(x)

A quantity per micron of length → per metre.
"""
permicron_to_perm(x::FT) where {FT} = x * FT(1e6)
permicron_to_perm(x::AbstractArray{FT}) where {FT <: AbstractFloat} = permicron_to_perm.(x)

"""Value `SHFOBS`/`LHFOBS` carry where they have no observation; not the declared `missing_value`."""
const ATLAS_UNDECLARED_SENTINEL = -99999.0

"""
    mask_undeclared_sentinel(x)

`NaN` where `x` is [`ATLAS_UNDECLARED_SENTINEL`](@ref), which nothing masks because it is not the
`missing_value` the file declares.
"""
mask_undeclared_sentinel(x::FT) where {FT} = x == FT(ATLAS_UNDECLARED_SENTINEL) ? FT(NaN) : x
mask_undeclared_sentinel(x::AbstractArray) = mask_undeclared_sentinel.(x)

"""
    ComposeFirst(f, g)

Callable applying `g` to the first argument only, then `f` to that result and any remaining
arguments: `ComposeFirst(f, g)(x, args...) == f(g(x), args...)`.
"""
struct ComposeFirst{F, G}
    f::F
    g::G
end
@inline (h::ComposeFirst)(x, args...) = h.f(h.g(x), args...)

# --- the spec ---------------------------------------------------------------------------------- #
"""
    AtlasVarSpec

How one Atlas variable is obtained.

  - `kind`        `:raw` (read from the file) or `:derived` (computed from other entries)
  - `transform`   `:raw` → `f(x, inputs...)`; `:derived` → `f(raw_inputs..., inputs...)`
  - `inputs`      other entries (already converted) the transform needs, e.g. `(:q_tot,)`
  - `raw_inputs`  file variables the transform needs UNCONVERTED
  - `raw_units`   asserted against the file's `units` attribute before converting
  - `units`       what the transform produces
  - `quantity`    the physical quantity, for grouping and for checking a consumer asked for the right thing
  - `provenance`  what settled the classification when the `units` attribute was not enough:
      * `:units_attr`         the attribute, and it is right
      * `:magnitude_verified` the attribute is absent or wrong; the values' magnitude settles it
      * `:long_name`          the attribute is ambiguous or wrong; `long_name` settles it
      * `:derived`            computed here rather than read
"""
struct AtlasVarSpec{F, I, R}
    kind::Symbol
    transform::F
    inputs::I
    raw_inputs::R
    raw_units::String
    units::String
    quantity::Symbol
    provenance::Symbol
end

"""Coordinates that are always available to a transform without being specs themselves."""
const ATLAS_COORDINATES = (:z, :time)

# --- variable groups --------------------------------------------------------------------------- #
# Generated from `ncdump -h` on the shipped archive, so the names are the file's own.

"""Microphysical process rates carrying a mass mixing ratio tendency, labelled `?/kg/day`."""
const ATLAS_MASS_PROCESS_RATES = (
    :EPRD, :EPRDG, :EPRDS, :EVPMG, :EVPMS, :MNUCCC,
    :MNUCCD, :MNUCCI, :MNUCCR, :PCC, :PCCN, :PGMLT,
    :PGRACS, :PGSACW, :PIACR, :PIACRS, :PITOSN, :PRA,
    :PRACG, :PRACI, :PRACIS, :PRACS, :PRAI, :PRC,
    :PRCI, :PRD, :PRDG, :PRDS, :PRE, :PSACR,
    :PSACWG, :PSACWI, :PSACWS, :PSMLT, :QHOMOC, :QHOMOR,
    :QMELTI, :QMULTG, :QMULTR, :QMULTRG, :QMULTS,
)

"""Microphysical process rates carrying a number tendency, labelled `?/kg/day`."""
const ATLAS_NUMBER_PROCESS_RATES = (
    :NGMLTG, :NGMLTR, :NGRACS, :NHOMOC, :NHOMOR, :NIACR,
    :NIACRS, :NMELTI, :NMULTG, :NMULTR, :NMULTRG, :NMULTS,
    :NNUCCC, :NNUCCD, :NNUCCI, :NNUCCR, :NPRA, :NPRACG,
    :NPRACS, :NPRAI, :NPRC, :NPRC1, :NPRCI, :NPSACWG,
    :NPSACWI, :NPSACWS, :NRAGG, :NSAGG, :NSCNG, :NSMLTR,
    :NSMLTS, :NSUBC, :NSUBG, :NSUBI, :NSUBR, :NSUBS,
)

"""Size-distribution slopes, labelled `?/kg/day` but inverse lengths (`long_name`: diameter = 1/(2λ))."""
const ATLAS_SLOPE_PARAMETERS = (:LAMC, :LAMG, :LAMI, :LAMR, :LAMS)

"""Dimensionless fractions and optical depths declaring empty units (`Ps` is the other empty one)."""
const ATLAS_EMPTY_UNIT_FRACTIONS = (
    :AREAPREC, :AREAPRTHR, :AUP, :CLD, :CLD245, :CLDCRM30,
    :CLDCRM40, :CLDCUMDN, :CLDCUMUP, :CLDHI, :CLDLOW, :CLDMID,
    :CLDSHD, :CLRAD00, :CLRAD05, :CLRAD10, :CLRAD15, :CLRAD20,
    :CLRADHI, :CLRADLO, :CLRADM05, :CLRADM10, :CLRADM15, :CLRADM20,
    :CLRADM25, :CLRADM30, :CLRADM35, :CORECL, :COREDNCL, :HYDRO,
    :ISCCPALB, :ISCCPHGH, :ISCCPLOW, :ISCCPMID, :ISCCPTAU, :ISCCPTOT,
    :MISRTOT, :MODISHGH, :MODISLOW, :MODISMID, :MODISTAU, :MODISTAUI,
    :MODISTAUL, :MODISTOT, :MODISTOTI, :MODISTOTL, :WSKEW,
)

"""Water fluxes reported as `W/m2` via the latent heat of vaporization."""
const ATLAS_LIQUID_WATER_FLUXES = (
    :QCFLUX, :QRFLXR, :QRFLXS, :QRSDFLX, :QTFLUX, :QTOFLXR,
    :QTOFLXS, :QTOSDFLX,
)

"""Water fluxes reported as `W/m2` via the latent heat of sublimation."""
const ATLAS_ICE_WATER_FLUXES = (
    :QGFLXR, :QGFLXS, :QGSDFLX, :QIFLUX, :QIFLXR, :QIFLXS,
    :QISDFLX, :QSFLXR, :QSFLXS, :QSSDFLX,
)

"""Genuine energy fluxes, already `W/m2`."""
const ATLAS_ENERGY_FLUXES = (
    :LHF, :LWDS, :LWNS, :LWNSC, :LWNT,
    :LWNTOA, :LWNTOAC, :RADLWDN, :RADLWUP, :RADSWDN, :RADSWUP,
    :SHF, :SOLIN, :SWDS, :SWNS, :SWNSC, :SWNT,
    :SWNTOA, :SWNTOAC, :TLFLUX, :TLFLUXS, :TVFLUX,
)

"""Observed fluxes carrying"""
const ATLAS_OBSERVED_ENERGY_FLUXES = (:LHFOBS, :SHFOBS)

"""Static energies stored divided by `cp_d`, so labelled `K`."""
const ATLAS_STATIC_ENERGIES = (:DSE, :DSECLD, :HFCLD, :HFCLDA, :MSE, :MSECLD, :SSE)

"""Temperatures and potential temperatures, already `K`."""
const ATLAS_TEMPERATURES = (
    :ISCCPTB, :ISCCPTBCLR, :SST, :SSTOBS, :TABS, :TABSOBS,
    :TACLD, :THETA, :THETAE, :THETAL, :THETAV, :TL,
    :TLCLD, :TVCLD, :TVCLDA,
)

"""Water mass mixing ratios in `g/kg`, converted with the local total water."""
const ATLAS_WATER_MASSES = (
    :QC, :QCI, :QCL, :QCOND, :QG, :QGCLD,
    :QI, :QICLD, :QN, :QNCLD, :QP, :QPCLD,
    :QPI, :QPL, :QR, :QRCLD, :QS, :QSAT,
    :QSCLD, :QT, :QTCLD, :QTO, :QTOCLD, :QV,
    :QVOBS,
)

"""Water mass mixing ratios mislabelled `kg/kg`; their magnitudes match their `g/kg` twins."""
const ATLAS_MISLABELLED_KGKG = (:QCCLD, :QVCLD)

"""Water mass mixing ratio tendencies in `g/kg/day`."""
const ATLAS_WATER_MASS_RATES = (
    :QGADV, :QGDIFF, :QGLSADV, :QGMPHY, :QGSED, :QHTEND,
    :QIADV, :QIDIFF, :QILSADV, :QIMPHY, :QISED, :QNUDGE,
    :QRADV, :QRDIFF, :QRLSADV, :QRMPHY, :QRSED, :QSADV,
    :QSDIFF, :QSLSADV, :QSMPHY, :QSSED, :QTEND, :QTOADV,
    :QTODIFF, :QTOLSADV, :QTOMPHY, :QTOSED, :QVTEND,
)

const ATLAS_NUMBER_CONCENTRATIONS = (
    :NCMN, :NG, :NGCLD, :NI, :NICLD, :NR,
    :NRCLD, :NRMN, :NS, :NSCLD,
)

const ATLAS_NUMBER_RATES = (
    :NGADV, :NGDIFF, :NGLSADV, :NGMPHY, :NGSED, :NIADV,
    :NIDIFF, :NILSADV, :NIMPHY, :NISED, :NRADV, :NRDIFF,
    :NRLSADV, :NRMPHY, :NRSED, :NSADV, :NSDIFF, :NSLSADV,
    :NSMPHY, :NSSED,
)

const ATLAS_NUMBER_FLUXES = (
    :NGFLXR, :NGFLXS, :NGSDFLX, :NIFLXR, :NIFLXS, :NISDFLX,
    :NRFLXR, :NRFLXS, :NRSDFLX, :NSFLXR, :NSFLXS, :NSSDFLX,
)

const ATLAS_TEMPERATURE_RATES = (
    :HLADV, :HLDIFF, :HLLAT, :HLRAD, :HLSTOR, :Q1C,
    :RADQR, :RADQRC, :RADQRLW, :RADQRS,
    :RADQRSW, :THTEND, :TNUDGE, :TTEND, :TVTEND,
)

"""Yanai apparent moisture sink and total-water storage: `K/day`, but `Q2 = -(L/cp) dq/dt`."""
const ATLAS_MOISTURE_RATES_IN_K = (:Q2, :QTSTOR)

"""Temperature tendencies spelling the same unit `K/d` rather than `K/day`."""
const ATLAS_RADQR_KDAY = (:RADQRCLW, :RADQRCSW)

const ATLAS_VELOCITY_RATES = (
    :UADV, :UDIFF, :ULSADVV, :UNUDGE, :URESID, :USTOR,
    :UTENDCOR, :VADV, :VDIFF, :VLSADVV, :VNUDGE, :VRESID,
    :VSTOR, :VTENDCOR,
)

const ATLAS_VELOCITIES = (
    :U, :UCLD, :UCLDA, :UMAX, :UOBS, :V,
    :VCLD, :VCLDA, :VOBS, :WCLD, :WCLDA, :WMAX,
    :WOBS,
)

const ATLAS_VELOCITY_MOMENT2 = (
    :TKE, :TKES, :U2, :UW, :UWCLD, :UWSB,
    :UWSBCLD, :V2, :VW, :VWCLD, :VWSB, :VWSBCLD,
    :W2,
)

const ATLAS_VELOCITY_MOMENT3 = (:W3, :WSTAR3)

const ATLAS_TKE_BUDGET = (
    :ADVTR, :ADVTRS, :BUOYA, :BUOYAS, :DIFTR, :DISSIP,
    :DISSIPS, :PRESSTR, :SHEAR, :SHEARS, :W2ADV, :W2BUOY,
    :W2DIFF, :W2PRES, :W2REDIS, :WUADV, :WUANIZ, :WUBUOY,
    :WUDIFF, :WUPRES, :WUSHEAR, :WVADV, :WVANIZ, :WVBUOY,
    :WVDIFF, :WVPRES, :WVSHEAR,
)

"""Condensed water paths in `g/m2`."""
const ATLAS_WATER_PATHS = (:CWP, :GWP, :IWP, :LWP, :MODISIWP, :MODISLWP, :RWP, :SWP)

const ATLAS_HEIGHTS_KM = (:MISRZTOP, :ZCB, :ZCBMIN, :ZCT, :ZCTMAX, :ZINV)
const ATLAS_AREAS_KM2 = (:ZCT2, :ZINV2)
const ATLAS_MASS_FLUXES = (:MC, :MCDNS, :MCDNU, :MCR, :MCRDNS, :MCRDNU, :MCRUP, :MCUP, :MFCLD)
const ATLAS_DIMENSIONLESS = (:QGFRAC, :QIFRAC, :QRFRAC, :QSFRAC, :TAUQC, :TAUQG, :TAUQI, :TAUQR, :TAUQS)
const ATLAS_SPECIFIC_ENERGIES = (:CAPE, :CAPEOBS, :CIN, :CINOBS)

"""Second moments of water mixing ratios in `g2/kg2`, taking `(1 - q_t)^2`."""
const ATLAS_WATER_MOMENT2 = (:QC2, :QI2, :QS2, :QT2)

const ATLAS_TEMP_MASS_FLUXES = (:MFTLCLD, :MFTLCLDA, :MFTVCLD, :MFTVCLDA)

"""Frozen-MSE mass fluxes: `K kg/m2/s`, but the static energy is stored divided by `cp_d`."""
const ATLAS_STATIC_ENERGY_MASS_FLUXES = (:MFHCLD, :MFHCLDA)
const ATLAS_TEMP_MOMENT2_RATES = (:T2ADVTR, :T2DIFTR, :T2DISSIP, :T2GRAD, :T2PREC)

"""Total-water variance budget terms, labelled `1/s` but carrying `(kg/kg)^2/s`."""
const ATLAS_WATER_VARIANCE_BUDGET = (:Q2ADVTR, :Q2DIFTR, :Q2DISSIP, :Q2GRAD, :Q2PREC)

"""Effective radii reported as `g/kg/micro`, a water content per micron of radius."""
const ATLAS_EFFECTIVE_RADII = (:QCOEFFR, :QGOEFFR, :QIOEFFR, :QROEFFR, :QSOEFFR)

const ATLAS_PRESSURES_MB = (:ISCCPPTOP, :MODISPTOP, :PRES, :p)
const ATLAS_MOMENTUM_FLUXES = (:RUWCLD, :RVWCLD, :RWWCLD)

"""Liquid-water-static-energy flux budget terms, labelled `m s-2 K`."""
const ATLAS_TEMP_FLUX_BUDGET = (:TWADV, :TWBUOY, :TWDIFF, :TWGRAD, :TWPREC, :TWPRES)

"""Total-water flux budget terms carrying a stray `K` in their label."""
const ATLAS_WATER_FLUX_BUDGET_K = (:QWADV, :QWDIFF, :QWGRAD)

"""Total-water flux budget terms, labelled `m s-2`."""
const ATLAS_WATER_FLUX_BUDGET = (:QWBUOY, :QWPREC, :QWPRES)

const ATLAS_WATER_VEL_COVAR = (:QCWCLD, :QIWCLD, :QTWCLD)
const ATLAS_TEMP_VEL_COVAR = (:TLWCLD, :TVWCLD)
const ATLAS_WATER_MASS_FLUX_G = (:MFQTCLD, :MFQTCLDA)
const ATLAS_ACCELERATIONS = (:UPGFCLD, :VPGFCLD, :WPGFCLD)
const ATLAS_DIFFUSIVITIES = (:TK, :TKH)
const ATLAS_PRECIP_RATES_DAY = (:PREC, :PRECIP)
const ATLAS_PRECIP_RATES_D = (:PRECMAX, :PRECMN)
const ATLAS_WATER_DEPTHS = (:PW, :PWOBS)
const ATLAS_MICRON_RADII = (:MODISREI, :MODISREL)

# --- the input registry ------------------------------------------------------------------------ #
"""
    atlas_input_var_specs()

Specs for the Atlas SAM *input* files, keyed by variable name.
"""
function atlas_input_var_specs(::Type{FT} = Float32) where {FT}
    raw(f, inputs, raw_units, units, quantity, provenance) =
        AtlasVarSpec(:raw, f, inputs, (), raw_units, units, quantity, provenance)

    FTB = Base.Fix1(broadcast, FT)
    return Dict{Symbol, AtlasVarSpec}(
        :lev => raw(FTB, (), "Pa", "Pa", :pressure, :units_attr),
        :Ps => raw(FTB, (), "Pa", "Pa", :pressure, :units_attr),
        :tsec => raw(FTB, (), "s", "s", :time, :units_attr),
        :T => raw(FTB, (), "K", "K", :temperature, :units_attr),
        :Tg => raw(FTB, (), "K", "K", :temperature, :units_attr),
        :q => raw(FTB ∘ wt_to_qt , (), "kg/kg", "kg/kg", :water_total_specific, :long_name), # Total water MASS MIXING RATIO. Atlas et al. (2020) define q_t as the sum of the mixing ratios of vapour and nonprecipitating condensate; the runs initialize cloud-free,
        :u => raw(FTB, (), "m/s", "m/s", :velocity, :units_attr),
        :v => raw(FTB, (), "m/s", "m/s", :velocity, :units_attr),
        :ug => raw(FTB, (), "m/s", "m/s", :velocity, :units_attr),
        :vg => raw(FTB, (), "m/s", "m/s", :velocity, :units_attr),
        :divT => raw(FTB, (), "K/s", "K/s", :temperature_rate, :units_attr),
        :divq => raw(FTB, (), "K/s", "kg/kg/s", :water_mass_rate, :long_name), # Declared "K/s", its long_name is the advection of a water vapour mixing ratio, so the quantity is a moisture rate and the label is simply wrong.
        :omega => raw(FTB, (), "Pa/s", "Pa/s", :pressure_rate, :units_attr),
        :Ptend => raw(FTB, (), "Pa/s", "Pa/s", :pressure_rate, :units_attr),
        :phis => raw(FTB, (), "m2/s2", "m2/s2", :geopotential, :units_attr),
        :o3mmr => raw(FTB, (), "kg/kg", "kg/kg", :ozone_mass, :units_attr),
        :lat => raw(FTB, (), "degree_north", "degree_north", :latitude, :units_attr),
        :lon => raw(FTB, (), "degree_east", "degree_east", :longitude, :units_attr),
    )
end

# --- the output registry ----------------------------------------------------------------------- #
"""
    atlas_output_var_specs(backend = DefaultThermodynamicsBackend(), FT = Float64; ρ_ice, q_small, sed_w_cap)

Specs for the Atlas LES *output* files, keyed by variable name, covering every variable the archive
carries plus the derived diagnostics.

`backend` supplies the latent heats and heat capacity the archive's reporting conventions imply, and
the saturation vapour pressures the humidity diagnostics need. `ρ_ice` sets the ice size
distribution, so it should be the density the archive was written with.
"""
function atlas_output_var_specs(
    backend::AbstractThermodynamicsBackend = DefaultThermodynamicsBackend(), # uses Atlas's constants by default
    ::Type{FT} = Float32; # default is native res (float 32)
    ρ_ice::FT = FT(ATLAS_ICE_APPARENT_DENSITY),
    q_small::FT = FT(ATLAS_Q_SMALL),
    sed_w_cap::@NamedTuple{cloud::FT, ice::FT, snow::FT, rain::FT} = map(FT, ATLAS_SED_W_CAP),
) where {FT}
    specs = Dict{Symbol, AtlasVarSpec}()

    FTB = Base.Fix1(broadcast, FT)

    function register!(names, kind, transform, inputs, raw_units, units, quantity, provenance)
        for name in names
            haskey(specs, name) && error("duplicate spec for Atlas variable $name")
            specs[name] = AtlasVarSpec(kind, transform, inputs, (), raw_units, units, quantity, provenance)
        end
        return nothing
    end
    raw!(names, f, raw_units, units, quantity, provenance; inputs = ()) =
        register!(names, :raw, f, inputs, raw_units, units, quantity, provenance)

    # The constants the LES encoded its fluxes with, not the backend's.
    L_v = L_v0(backend, FT) # latent_heat_vapor(thermo_params, T) is better, but SAM is fixed Lv
    L_s = L_s0(backend, FT) # latent_heat_sublim(thermo_params, T) is better, but SAM is fixed Ls
    cp = cp_d(backend, FT)

    # ---- coordinates ------------------------------------------------------------------------- #
    raw!((:z,), identity, "m", "m", :height, :units_attr)
    raw!((:time,), day_to_s, "day", "s", :time, :units_attr)

    # ---- pressure ---------------------------------------------------------------------------- #
    raw!(ATLAS_PRESSURES_MB, hPa_to_Pa, "mb", "Pa", :pressure, :units_attr)
    # Declares no units at all; it reads ~984 on RF09, which is hPa rather than Pa.
    raw!((:Ps,), FTB ∘ hPa_to_Pa, "", "Pa", :pressure, :magnitude_verified)

    # ---- density ----------------------------------------------------------------------------- #
    raw!((:RHO,), FTB, "kg/m3", "kg/m3", :density, :units_attr)

    # ---- temperature ------------------------------------------------------------------------- #
    raw!(ATLAS_TEMPERATURES, FTB, "K", "K", :temperature, :units_attr)
    raw!(ATLAS_STATIC_ENERGIES, FTB ∘ K_to_Jkg, "K", "J/kg", :specific_energy, :long_name; inputs = (cp,))
    raw!(ATLAS_TEMPERATURE_RATES, FTB ∘ perday_to_persec, "K/day", "K/s", :temperature_rate, :units_attr)
    raw!(ATLAS_RADQR_KDAY, FTB ∘ perday_to_persec, "K/d", "K/s", :temperature_rate, :units_attr)
    # Yanai apparent moisture sink and total-water storage: reported as K/day via Q2 = -(L/cp) dq/dt.
    raw!(
        ATLAS_MOISTURE_RATES_IN_K, FTB ∘ ComposeFirst(wc_to_qc, ComposeFirst(Base.Fix1(*, cp / L_v), perday_to_persec)),
        "K/day", "kg/kg/s", :water_mass_rate, :long_name; inputs = (:q_tot,),
    )
    raw!((:TL2,), FTB, "K2", "K2", :temperature_moment2, :units_attr)
    raw!(ATLAS_TEMP_MOMENT2_RATES, FTB, "K2/s", "K2/s", :temperature_moment2_rate, :units_attr)
    raw!(ATLAS_SPECIFIC_ENERGIES, FTB, "J/kg", "J/kg", :specific_energy, :units_attr)

    # ---- water masses ------------------------------------------------------------------------ #
    # A mixing ratio becomes a specific content by the local total water, so every one of these
    # takes `q_tot` (itself derived below from QT + QP).
    raw!(
        ATLAS_WATER_MASSES, FTB ∘ ComposeFirst(wc_to_qc, g_to_kg), "g/kg", "kg/kg",
        :water_mass, :units_attr; inputs = (:q_tot,),
    )
    raw!(
        ATLAS_MISLABELLED_KGKG, FTB ∘ ComposeFirst(wc_to_qc, g_to_kg), "kg/kg", "kg/kg",
        :water_mass, :magnitude_verified; inputs = (:q_tot,),
    )
    raw!(
        ATLAS_WATER_MASS_RATES, FTB ∘ ComposeFirst(wc_to_qc, ComposeFirst(g_to_kg, perday_to_persec)),
        "g/kg/day", "kg/kg/s", :water_mass_rate, :units_attr; inputs = (:q_tot,),
    )
    raw!(
        ATLAS_WATER_MOMENT2, FTB ∘ ComposeFirst(wc2_to_qc2, g2_to_kg2), "g2/kg2", "kg2/kg2",
        :water_moment2, :units_attr; inputs = (:q_tot,),
    )
    raw!(ATLAS_WATER_PATHS, FTB ∘ g_to_kg, "g/m2", "kg/m2", :water_path, :units_attr)
    raw!((:LWP2,), FTB ∘ g2_to_kg2, "(g/m2)^2", "kg2/m4", :water_path_moment2, :units_attr)
    raw!((:QTO2,), FTB ∘ ComposeFirst(wc2_to_qc2, g2_to_kg2), "(g/kg)^2", "kg2/kg2", :water_moment2, :units_attr; inputs = (:q_tot,))
    # Precipitable water is a depth of liquid; 1 mm over 1 m2 is ρ_liq/1000 kg. Atlas ρ_liq = 1000
    raw!(ATLAS_WATER_DEPTHS, FTB ∘ mm_to_m ∘ Base.Fix2((/), 1000.0), "mm", "kg/m2", :water_path, :long_name)

    # ---- microphysical process rates --------------------------------------------------------- #
    # All three groups declare `?/kg/day`;
    raw!(
        ATLAS_MASS_PROCESS_RATES, FTB ∘ ComposeFirst(wc_to_qc, perday_to_persec),
        "?/kg/day", "kg/kg/s", :water_mass_rate, :long_name; inputs = (:q_tot,),
    )
    raw!(ATLAS_NUMBER_PROCESS_RATES, FTB ∘ perday_to_persec, "?/kg/day", "1/kg/s", :number_mass_rate, :long_name)
    raw!(ATLAS_SLOPE_PARAMETERS, FTB ∘ identity, "?/kg/day", "1/m", :inverse_length, :long_name)

    # ---- number concentrations, rates and fluxes --------------------------------------------- #
    raw!(ATLAS_NUMBER_CONCENTRATIONS, FTB ∘ percm3_to_perm3, "#/cm3", "1/m3", :number_concentration, :units_attr)
    raw!(
        ATLAS_NUMBER_RATES, FTB ∘ ComposeFirst(percm3_to_perm3, perday_to_persec),
        "#/cm3/day", "1/m3/s", :number_rate, :units_attr,
    )
    raw!(ATLAS_NUMBER_FLUXES, FTB, "#/m2/s", "1/m2/s", :number_flux, :units_attr)

    # ---- fluxes ------------------------------------------------------------------------------ #
    raw!(ATLAS_ENERGY_FLUXES, FTB, "W/m2", "W/m2", :energy_flux, :units_attr)
    # Both carry -99999 where there is no observation, which is not the declared missing_value.
    raw!((:LHFOBS,), FTB ∘ mask_undeclared_sentinel, "W/m2", "W/m2", :energy_flux, :magnitude_verified)
    # The observed sensible heat flux declares its own name as its units.
    raw!((:SHFOBS,), FTB ∘ mask_undeclared_sentinel, "SHFOBS", "W/m2", :energy_flux, :long_name)
    # Water fluxes reported as energy fluxes; the phase sets which latent heat was used.
    raw!(
        ATLAS_LIQUID_WATER_FLUXES, FTB ∘ wm2_to_kgm2s, "W/m2", "kg/m2/s",
        :water_mass_flux, :long_name; inputs = (L_v,),
    )
    raw!(
        ATLAS_ICE_WATER_FLUXES, FTB ∘ wm2_to_kgm2s, "W/m2", "kg/m2/s",
        :water_mass_flux, :long_name; inputs = (L_s,),
    )
    raw!(ATLAS_MASS_FLUXES, FTB, "kg/m2/s", "kg/m2/s", :mass_flux, :units_attr)
    raw!(ATLAS_TEMP_MASS_FLUXES, FTB, "K kg/m2/s", "K kg/m2/s", :temperature_mass_flux, :units_attr)
    raw!(
        ATLAS_STATIC_ENERGY_MASS_FLUXES, FTB ∘ K_to_Jkg, "K kg/m2/s", "W/m2",
        :static_energy_mass_flux, :long_name; inputs = (cp,),
    )
    raw!(ATLAS_WATER_MASS_FLUX_G, FTB ∘ ComposeFirst(wc_to_qc, g_to_kg), "g/m2/s", "kg/m2/s", :water_mass_flux, :units_attr; inputs = (:q_tot,))
    raw!(ATLAS_MOMENTUM_FLUXES, identity, "kg/m/s2", "kg/m/s2", :momentum_flux, :units_attr)

    # ---- velocities and their moments -------------------------------------------------------- #
    raw!(ATLAS_VELOCITIES, FTB, "m/s", "m/s", :velocity, :units_attr)
    raw!(ATLAS_VELOCITY_RATES, FTB ∘ perday_to_persec, "m/s/day", "m/s2", :acceleration, :units_attr)
    raw!(ATLAS_ACCELERATIONS, FTB, "m/s2", "m/s2", :acceleration, :units_attr)
    raw!(ATLAS_VELOCITY_MOMENT2, FTB, "m2/s2", "m2/s2", :velocity_moment2, :units_attr)
    raw!(ATLAS_VELOCITY_MOMENT3, FTB, "m3/s3", "m3/s3", :velocity_moment3, :units_attr)
    raw!(ATLAS_TKE_BUDGET, FTB, "m2/s3", "m2/s3", :tke_budget, :units_attr)
    raw!(ATLAS_DIFFUSIVITIES, FTB, "m2/s", "m2/s", :diffusivity, :units_attr)

    # ---- covariance and flux budgets --------------------------------------------------------- #
    raw!(ATLAS_TEMP_VEL_COVAR, FTB, "Km/s", "K m/s", :temperature_velocity_covar, :units_attr)
    raw!(
        ATLAS_WATER_VEL_COVAR, FTB ∘ ComposeFirst(wc_to_qc, g_to_kg), "g/kg m/s", "kg/kg m/s",
        :water_velocity_covar, :units_attr; inputs = (:q_tot,),
    )
    # Covariance of HL [K] and QT, already K*kg/kg, so it takes one water factor and no g_to_kg.
    raw!((:TQ,), FTB ∘ wc_to_qc, "K2", "K kg/kg", :temperature_water_covar, :magnitude_verified; inputs = (:q_tot,))
    raw!(ATLAS_TEMP_FLUX_BUDGET, FTB, "m s-2 K", "K m/s2", :temperature_flux_budget, :units_attr)
    raw!(
        ATLAS_WATER_FLUX_BUDGET_K, FTB ∘ wc_to_qc, "m s-2 K", "m/s2 kg/kg",
        :water_flux_budget, :magnitude_verified; inputs = (:q_tot,),
    )
    raw!(
        ATLAS_WATER_FLUX_BUDGET, FTB ∘ wc_to_qc, "m s-2", "m/s2 kg/kg",
        :water_flux_budget, :magnitude_verified; inputs = (:q_tot,),
    )
    raw!(
        ATLAS_WATER_VARIANCE_BUDGET, FTB ∘ ComposeFirst(wc2_to_qc2, g2_to_kg2), "1/s", "kg2/kg2/s",
        :water_moment2_rate, :magnitude_verified; inputs = (:q_tot,),
    )

    # ---- lengths, areas, fractions ----------------------------------------------------------- #
    raw!(ATLAS_HEIGHTS_KM, FTB ∘ km_to_m, "km", "m", :height, :units_attr)
    raw!(ATLAS_AREAS_KM2, FTB ∘ km2_to_m2, "km2", "m2", :height_moment2, :units_attr)
    # Labelled `km`, but a variance of a height cannot be a length; its siblings say km2.
    raw!((:ZCB2,), FTB ∘ km2_to_m2, "km", "m2", :height_moment2, :long_name)
    raw!(ATLAS_MICRON_RADII, FTB ∘ micron_to_m, "mkm", "m", :length, :magnitude_verified)
    raw!(ATLAS_EMPTY_UNIT_FRACTIONS, FTB, "", "1", :dimensionless, :long_name)
    raw!(ATLAS_DIMENSIONLESS, FTB, "1", "1", :dimensionless, :units_attr)
    raw!((:RELH,), FTB ∘ percent_to_fraction, "per cent", "1", :relative_humidity, :units_attr)

    # ---- precipitation ----------------------------------------------------------------------- #
    # 1 mm/day of liquid over 1 m2 is 1 kg/m2/day.
    raw!(ATLAS_PRECIP_RATES_DAY, FTB ∘ perday_to_persec, "mm/day", "kg/m2/s", :water_mass_flux, :long_name)
    raw!(ATLAS_PRECIP_RATES_D, FTB ∘ perday_to_persec, "mm/d", "kg/m2/s", :water_mass_flux, :long_name)
    raw!(
        (:PREC2,), FTB ∘ perday_to_persec ∘ perday_to_persec, "(mm/d)^2", "kg2/m4/s2",
        :water_mass_flux_moment2, :long_name,
    )

    # ---- effective radii --------------------------------------------------------------------- #
    raw!(
        ATLAS_EFFECTIVE_RADII, FTB ∘ ComposeFirst(wc_to_qc, permicron_to_perm ∘ g_to_kg),
        "g/kg/micro", "kg/kg/m", :water_mass_per_length, :units_attr; inputs = (:q_tot,),
    )

    add_atlas_derived_specs!(specs, backend, FT; ρ_ice, q_small, sed_w_cap)
    return specs
end

"""
    add_atlas_derived_specs!(specs, backend, FT; ρ_ice, q_small, sed_w_cap)

Add the derived entries to `specs`, in place. Their `inputs` are entries of `specs`, so they arrive
already converted to the units those entries declare.
"""
function add_atlas_derived_specs!(
    specs::AbstractDict{Symbol, AtlasVarSpec},
    backend::AbstractThermodynamicsBackend,
    ::Type{FT} = Float32;
    ρ_ice::FT,
    q_small::FT,
    sed_w_cap::@NamedTuple{cloud::FT, ice::FT, snow::FT, rain::FT},
) where {FT}
    function derived!(
        name::Symbol, transform, inputs::Tuple, units::AbstractString, quantity::Symbol;
        raw_inputs::Tuple = (),
    )
        haskey(specs, name) && error("duplicate spec for Atlas variable $name")
        specs[name] = AtlasVarSpec(:derived, transform, inputs, raw_inputs, "", units, quantity, :derived)
        return nothing
    end

    # ---- total water ------------------------------------------------------------------------- #
    # `QT` is vapour + cloud liquid + cloud ice; rain and snow sit in `QP`, and SAM's buoyancy
    # carries `qpl+qpi` in the loading term, so both are part of the parcel mass
    derived!(
        :q_tot, (qt, qp) -> wt_to_qt(g_to_kg(qt) .+ g_to_kg(qp)), (), "kg/kg", :water_total_specific;
        raw_inputs = (:QT, :QP),
    )

    # ---- standard deviations of the water moments -------------------------------------------- #
    for (std_name, var_name) in ((:QT2_STD, :QT2), (:QC2_STD, :QC2), (:QI2_STD, :QI2), (:QS2_STD, :QS2))
        derived!(std_name, x -> sqrt.(x), (var_name,), "kg/kg", :water_mass)
    end

    # ---- aliases ------------------------------------------------------------------------------ #
    # `QTO` is vapour + cloud liquid and vapour does not sediment, so the sedimentation of QTO is
    # that of cloud liquid alone.
    derived!(:QCSED, identity, (:QTOSED,), "kg/kg/s", :water_mass_rate)
    derived!(:QCSDFLX, identity, (:QTOSDFLX,), "kg/m2/s", :water_mass_flux)

    # ---- net deposition rates ----------------------------------------------------------------- #
    derived!(:PRDT, (+), (:PRD, :EPRD), "kg/kg/s", :water_mass_rate)
    derived!(:PRDST, (+), (:PRDS, :EPRDS), "kg/kg/s", :water_mass_rate)

    # ---- relative humidity -------------------------------------------------------------------- #
    for (name, q_sat) in ((:RHL, q_vap_saturation_liq), (:RHI, q_vap_saturation_ice))
        derived!(
            name, (qv, T, p) -> qv ./ q_sat.(Ref(backend), T, p),
            (:QV, :TABS, :PRES), "1", :relative_humidity,
        )
        # Bracketing values one total-water standard deviation either side of the vapour.
        for (suffix, sgn) in ((:_low, -one(FT)), (:_high, one(FT)))
            derived!(
                Symbol(name, suffix), (qv, T, p, σ) -> (qv .+ sgn .* σ) ./ q_sat.(Ref(backend), T, p),
                (:QV, :TABS, :PRES, :QT2_STD), "1", :relative_humidity,
            )
        end
    end

    # ---- sedimentation velocities -------------------------------------------------------------- #
    for (name, flux, species, cap) in (
        (:WI, :QISDFLX, :QI, sed_w_cap.ice),
        (:WS, :QSSDFLX, :QS, sed_w_cap.snow),
        (:WR, :QRSDFLX, :QR, sed_w_cap.rain),
        (:WC, :QCSDFLX, :QC, sed_w_cap.cloud),
    )
        derived!(
            name, (F, q, ρ) -> sedimentation_velocity.(F, q, ρ, cap),
            (flux, species, :RHO), "m/s", :velocity,
        )
    end

    # ---- ice number ---------------------------------------------------------------------------- #
    derived!(:NI_ALL, (a, b) -> a .+ b, (:NI, :NS), "1/m3", :number_concentration)
    derived!(:NINP, T -> cooper_ice_nucleus_concentration.(T), (:TABS,), "1/m3", :number_concentration)

    # ---- ice particle radius and relaxation timescales ------------------------------------------ #
    derived!(:RI, (q, N, ρ) -> ice_mean_radius.(q, N, ρ, ρ_ice), (:QI, :NI, :RHO), "m", :length)

    # `_full` is the whole ice distribution; the bare name is the part below the 125 µm threshold,
    # which is the part Atlas reports a deposition rate for. A smaller population takes up less
    # vapour, hence `τ_below = τ_full / weight >= τ_full`.
    derived!(
        :τ_ice_theory_full, (q, N, ρ, T, p) -> ice_theory_timescale.(q, N, ρ, T, p, ρ_ice),
        (:QI, :NI, :RHO, :TABS, :PRES), "s", :timescale,
    )
    derived!(
        :τ_ice_theory, (τ, q, N, ρ) -> τ ./ ice_process_threshold_weight.(q, N, ρ, ρ_ice),
        (:τ_ice_theory_full, :QI, :NI, :RHO), "s", :timescale,
    )

    derived!(
        :τ_ice_full,
        function (qv, T, p, ρ, prdt, q_sno, q, N)
            excess = qv .- q_vap_saturation_ice.(Ref(backend), T, p)
            below = ice_process_threshold_weight.(q, N, ρ, ρ_ice)
            rate = ifelse.(q_sno .>= q_small, prdt ./ below, prdt)
            return excess ./ ifelse.(rate .== 0, FT(NaN), rate)
        end,
        (:QV, :TABS, :PRES, :RHO, :PRDT, :QS, :QI, :NI), "s", :timescale,
    )
    derived!(
        :τ_ice, (τ, q, N, ρ) -> τ ./ ice_process_threshold_weight.(q, N, ρ, ρ_ice),
        (:τ_ice_full, :QI, :NI, :RHO), "s", :timescale,
    )

    return specs
end

"""
    atlas_var_specs(source, backend = DefaultThermodynamicsBackend(), FT = Float64; kwargs...)

The registry for a regrid `source`: [`AtlasInput`](@ref) or [`LESOutput`](@ref).
"""
atlas_var_specs(::AtlasInput, ::Type{FT} = Float32) where {FT <: AbstractFloat} = atlas_input_var_specs(FT)
atlas_var_specs(::LESOutput,
    backend::DefaultThermodynamicsBackend = DefaultThermodynamicsBackend(), 
    ::Type{FT} = Float32; 
kwargs...) where {FT <: AbstractFloat} = atlas_output_var_specs(backend, FT; kwargs...)

# --- dependency order --------------------------------------------------------------------------- #
"""
    atlas_processing_order(specs; seeds = ATLAS_COORDINATES)

Entries ordered so every `inputs` dependency precedes its consumer, by Kahn's algorithm. Ties break
by name, so the order and any failure are reproducible run to run. Errors on a dependency cycle,
naming the entries left unordered.
"""
function atlas_processing_order(specs::AbstractDict{Symbol, AtlasVarSpec}; seeds = ATLAS_COORDINATES)
    deps(spec) = Tuple(d for d in spec.inputs if d isa Symbol)
    for (name, spec) in specs, dep in deps(spec)
        (haskey(specs, dep) || dep in seeds) ||
            error("input $dep of $name is neither a spec nor a seeded coordinate")
    end

    indegree = Dict{Symbol, Int}(name => count(d -> !(d in seeds), deps(spec)) for (name, spec) in specs)
    dependents = Dict{Symbol, Vector{Symbol}}(name => Symbol[] for name in keys(specs))
    for (name, spec) in specs, dep in deps(spec)
        dep in seeds || push!(dependents[dep], name)
    end

    ready = sort!([name for (name, d) in indegree if d == 0])
    order = Symbol[]
    while !isempty(ready)
        name = popfirst!(ready)
        push!(order, name)
        for consumer in dependents[name]
            indegree[consumer] -= 1
            if indegree[consumer] == 0
                insert!(ready, searchsortedfirst(ready, consumer), consumer)
            end
        end
    end
    length(order) == length(specs) || error(
        "cyclic dependency among Atlas specs; unordered: $(sort!(setdiff(collect(keys(specs)), order)))",
    )
    return order
end

# --- reading through the registry ---------------------------------------------------------------- #
"""
    declared_units(var)

The `units` attribute of a NetCDF variable, or `nothing` when it carries none. An empty attribute is
reported as `""`, which is a declaration of nothing rather than an absence, and specs distinguish the
two.
"""
declared_units(var) = hasproperty(var, :attrib) ? get(var.attrib, "units", nothing) : nothing


"""
    convert_to_SI_units(data, name, source; specs = atlas_var_specs(source), inputs = (;))

Read `data[name]` and convert it to the units its spec declares
"""
function convert_to_SI_units(
    data,
    name::Union{Symbol, AbstractString},
    source::AbstractRegridSource,
    ::Type{FT} = Base.nonmissingtype(eltype(data[name]))
    ;
    backend = DefaultThermodynamicsBackend(),
    specs = (source isa AtlasInput) ? atlas_input_var_specs(FT) : atlas_output_var_specs(backend, FT),
    inputs::NamedTuple = (;),
) where {FT}
    key = Symbol(name)
    spec = specs[key]

    var = data[name]
    FTR = (FT <: Nothing) ? Base.nonmissingtype(eltype(var)) : FT
    var = (eltype(var) isa FTR) ? var : convert(Array{FTR}, var)

    args = map(spec.inputs) do input
        input isa Symbol || return input
        return inputs[input]
    end
    return spec.transform(_materialize(var), args...)
end

"""
    read_atlas_variable(data, name, source; specs = atlas_var_specs(source), cache = Dict{Symbol, Any}())

`data[name]` in the units its spec declares. A `:raw` entry is read and converted; a `:derived` entry
is computed from its `inputs`, resolved recursively so they arrive converted, and its `raw_inputs`,
read from the file unconverted.

`cache` memoizes across the call tree and may be reused across calls on the same `data`, so a shared
dependency like `q_tot` is read and converted once.
"""
function read_atlas_variable(
    data,
    name::Union{Symbol, AbstractString},
    source::AbstractRegridSource,
    backend::DefaultThermodynamicsBackend = DefaultThermodynamicsBackend(),
    ::Type{FT} = Float32;
    specs = (source isa AtlasInput) ? atlas_input_var_specs(FT) : atlas_output_var_specs(backend, FT),
    cache::Dict{Symbol, Any} = Dict{Symbol, Any}(),
) where {FT}
    key = Symbol(name)
    haskey(cache, key) && return cache[key]
    spec = specs[key]

    FTR = (FT <: Nothing) ? Base.nonmissingtype(eltype(data[String(name)])) : FT

    resolved = (;
        (n => read_atlas_variable(data, n, source; specs, cache) for n in spec.inputs if n isa Symbol)...
    )
    value = if spec.kind === :raw
        convert_to_SI_units(data, key, source, FTR; specs, inputs = resolved)
    else
        raw_args = map(n -> _materialize(data[String(n)], FTR), spec.raw_inputs)
        args = map(n -> n isa Symbol ? resolved[n] : n, spec.inputs)
        spec.transform(raw_args..., args...)
    end
    cache[key] = value
    return value
end