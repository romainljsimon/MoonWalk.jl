module MoonWalk

using Random, Distributions, Plots, LinearAlgebra, StaticArrays
using LaTeXStrings, NaNMath, DSP
using Plots.PlotMeasures

export RotationParameters, InfinityInteger, solve_sde

struct InfinityInteger <: Integer end

mutable struct RotationParameters{I<:Integer}
    dt::Float64
    T::Float64
    walkers::Int
    order::I
    simulation::String
    H::Float64
    rate::Float64
    cage_time::Float64
    escape_time::Float64
end

# Constructors
function RotationParameters(dt::Float64, T::Float64, simulation::String, walkers::Int, order::Integer, H::Float64, rate::Float64, cage_time::Float64, escape_time::Float64)
    RotationParameters(dt, T, walkers, order, simulation, H, rate, cage_time, escape_time)
end

function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer)
    RotationParameters(dt, T, walkers, order, "Brownian", 0.0, 0.0, 0.0, Inf)
end

function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer, H::Float64)
    RotationParameters(dt, T, walkers, order, "Cage", H, 0.0, Inf, 0.0)
end

function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer, H::Float64, rate::Float64)
    RotationParameters(dt, T, walkers, order, "Escape", H, rate, 0.0, 0.0)
end

# Utility Functions
include("utils.jl")

end # module