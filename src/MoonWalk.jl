module MoonWalk

using Random, Distributions, Plots, LinearAlgebra, StaticArrays, ProgressMeter, Printf
using NaNMath

export BrownianParameters, CageParameters, ParetoParameters
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
    amplitude::Float64
    α::Float16
end


# Utility Functions
include("angle-definition.jl")
include("io.jl")
include("math_utils.jl")
include("simulation.jl")

end # module
