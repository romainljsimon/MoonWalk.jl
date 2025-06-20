module MoonWalk

using Random, Distributions, Plots, LinearAlgebra, StaticArrays
using NaNMath

export RotationParameters, InfinityInteger, solve_sde

"""
    InfinityInteger

A type used to represent an infinite integer, mainly for specifying infinite order in BCH expansions.
"""
struct InfinityInteger <: Integer end

"""
    RotationParameters

Structure holding all parameters required for simulating rotational dynamics.

# Fields
- `dt::Float64`: Time step of the simulation.
- `T::Float64`: Total simulation time.
- `walkers::Int`: Number of walkers (simulated trajectories).
- `order::I`: Order of the BCH expansion (can be an integer or `InfinityInteger`).
- `simulation::String`: Type of simulation ("Brownian", "Cage", "Escape").
- `H::Float64`: Cage height (for "Cage" or "Escape" simulations).
- `rate::Float64`: Escape rate (for "Escape").
- `cage_time::Float64`: Time spent in the cage (for "Cage").
- `escape_time::Float64`: Escape time (for "Escape").
"""
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

"""
    RotationParameters(dt, T, simulation, walkers, order, H, rate, cage_time, escape_time)

General constructor for the `RotationParameters` structure.
"""
function RotationParameters(dt::Float64, T::Float64, simulation::String, walkers::Int, order::Integer, H::Float64, rate::Float64, cage_time::Float64, escape_time::Float64)
    RotationParameters(dt, T, walkers, order, simulation, H, rate, cage_time, escape_time)
end

"""
    RotationParameters(dt, T, walkers, order)

Constructor for a simple Brownian simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer)
    RotationParameters(dt, T, walkers, order, "Brownian", 0.0, 0.0, 0.0, Inf)
end

"""
    RotationParameters(dt, T, walkers, order, H)

Constructor for a "Cage" simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer, H::Float64)
    RotationParameters(dt, T, walkers, order, "Cage", H, 0.0, Inf, 0.0)
end

"""
    RotationParameters(dt, T, walkers, order, H, rate)

Constructor for an "Escape" simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int, order::Integer, H::Float64, rate::Float64)
    RotationParameters(dt, T, walkers, order, "Escape", H, rate, 0.0, 0.0)
end

# Utility Functions
include("math_utils.jl")
include("simulation.jl")

end # module