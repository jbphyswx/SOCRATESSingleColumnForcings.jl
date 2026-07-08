#=

=#


# Singleton BC types
abstract type AbstractBoundaryCondition end
struct ExtrapolateBoundaryCondition <: AbstractBoundaryCondition end
struct ErrorBoundaryCondition <: AbstractBoundaryCondition end
struct NearestBoundaryCondition <: AbstractBoundaryCondition end


const ValidBoundaryConditions = Union{ExtrapolateBoundaryCondition, ErrorBoundaryCondition, NearestBoundaryCondition}


bc_string(::ExtrapolateBoundaryCondition) = "extrapolate"
bc_string(::ErrorBoundaryCondition) = "error"
bc_string(::NearestBoundaryCondition) = "nearest"


bc_symbol(::ExtrapolateBoundaryCondition) = :extrapolate
bc_symbol(::ErrorBoundaryCondition) = :error
bc_symbol(::NearestBoundaryCondition) = :nearest


create_boundary_condition(x::AbstractBoundaryCondition) = x # pass through
create_boundary_condition(s::String) = create_boundary_condition(Symbol(lowercase(strip(s))))
create_boundary_condition(s::Symbol) = create_boundary_condition(Val(s))

# valid BCs
create_boundary_condition(::Val{:extrapolate}) = ExtrapolateBoundaryCondition()
create_boundary_condition(::Val{:error}) = ErrorBoundaryCondition()
create_boundary_condition(::Val{:nearest}) = NearestBoundaryCondition()

# public constructor alias: accepts a BC instance, Symbol, or String and funnels to
# `create_boundary_condition`; an unknown Symbol/String hits the `Val{s}` error below.
create_bc(x) = create_boundary_condition(x)

# fallback (dynamic, no extra functions)


create_bc(::Val{s}) where {s} = error(
    "Invalid boundary condition symbol: $s. Valid options are " *
    repr((T -> bc_string(T())).(Base.uniontypes(ValidBoundaryConditions))),
)