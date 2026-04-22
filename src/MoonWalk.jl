module MoonWalk

using Random, Distributions, LinearAlgebra, StaticArrays

export BrownianParameters, CageParameters, ParetoParameters, CageEscapeParameters
export simulation
export ExactRotation, IntegralDefinition, UnboundedDefinition


"""
    RotationParameters

Structure holding all parameters required for simulating rotational dynamics.
"""

abstract type RotationParameters end

struct BrownianParameters <: RotationParameters
    T::Float64
    amplitude::Float64
end

struct CageParameters <: RotationParameters
    T::Float64
    amplitude::Float64
    cage_size::Float64
end

struct ParetoParameters <: RotationParameters
    T::Float64
    amplitude_small::Float64
    amplitude_large::Float64
    α::Float16
    cage_size::Float64
    τ::Float64
end

struct CageEscapeParameters <: RotationParameters
    T::Float64
    amplitude_small::Float64
    amplitude_large::Float64
    rate::Float64
    cage_size::Float64
end


# Utility Functions
include("angle-definition.jl")
include("io.jl")
include("math_utils.jl")
include("simulation.jl")

end # module
