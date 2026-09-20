#=

=#


# Singleton BC types
"""Supertype for out-of-range evaluation policies on 1D interpolants."""
abstract type AbstractBoundaryCondition end

"""Linear extrapolation outside the node range."""
struct ExtrapolateBoundaryCondition <: AbstractBoundaryCondition end

"""Error when evaluating outside the node range."""
struct ErrorBoundaryCondition <: AbstractBoundaryCondition end

"""Return the nearest endpoint value outside the node range."""
struct NearestBoundaryCondition <: AbstractBoundaryCondition end

"""Return a constant value outside the node range."""
struct ConstantBoundaryCondition{T} <: AbstractBoundaryCondition
    value::T
end


const ValidBoundaryConditions = Union{ExtrapolateBoundaryCondition, ErrorBoundaryCondition, NearestBoundaryCondition, ConstantBoundaryCondition}


bc_string(::ExtrapolateBoundaryCondition) = "extrapolate"
bc_string(::ErrorBoundaryCondition) = "error"
bc_string(::NearestBoundaryCondition) = "nearest"
bc_string(x::ConstantBoundaryCondition; full::Bool=false) = full ? "constant($(x.value))" : "constant"


bc_symbol(::ExtrapolateBoundaryCondition) = :extrapolate
bc_symbol(::ErrorBoundaryCondition) = :error
bc_symbol(::NearestBoundaryCondition) = :nearest
bc_symbol(x::ConstantBoundaryCondition; full::Bool=false) = Symbol(bc_string(x; full=full))


create_boundary_condition(x::AbstractBoundaryCondition) = x # pass through
create_boundary_condition(s::String) = create_boundary_condition(Symbol(lowercase(strip(s))))
create_boundary_condition(s::Symbol) = create_boundary_condition(Val(s))

# valid BCs
create_boundary_condition(::Val{:extrapolate}) = ExtrapolateBoundaryCondition()
create_boundary_condition(::Val{:error}) = ErrorBoundaryCondition()
create_boundary_condition(::Val{:nearest}) = NearestBoundaryCondition()
create_boundary_condition(::Val{:constant}, value::FT) where {FT} = ConstantBoundaryCondition(value)
function create_boundary_condition(::Val{s}) where {s}
    str = String(s)
    if startswith(str, "constant(")
        return ConstantBoundaryCondition(Meta.parse(str[(length("constant(") + 1):(end - 1)]))
    end
    error("Invalid boundary condition: $s. Valid options are :extrapolate, :error, :nearest, \"constant(<value>)\"")
end

# public constructor alias: accepts a BC instance, Symbol, or String and funnels to
# `create_boundary_condition`.
"""
    create_bc(x)

Create a boundary condition from a singleton instance, `Symbol`, or `String`
(`:extrapolate`, `:error`, `:nearest`).
"""
create_bc(x) = create_boundary_condition(x)