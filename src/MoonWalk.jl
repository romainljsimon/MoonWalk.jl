module MoonWalk

using Random, Distributions, Plots, LinearAlgebra, StaticArrays, ProgressMeter, Printf, JLD2
using NaNMath

export RotationParameters, InfinityInteger, simulation, load_timestep
export get_time_trajectory, create_log_scheduler, load_trajectory_walker
export analyze_trajectory, ThetaMethod, IntegralMethod, ThreshThetaMethod, EulerMethod, UnboundedThetaMethod
export order0, order1, order2, order3
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
- `simulation::String`: Type of simulation ("Brownian", "Cage", "Escape").
- `H::Float64`: Cage height (for "Cage" or "Escape" simulations).
- `rate::Float64`: Escape rate (for "Escape").
- `tᵪ::Float64`: Time spent in the cage (for "Cage").
- `tₑ::Float64`: Escape time (for "Escape").
"""
mutable struct RotationParameters{D}
    dt::Float64
    T::Float64
    walkers::Int
    simulation::String
    Dᵣ::Float64
    H::Float64
    rate::Float64
    tᵪ::MVector{D, Float64}
    tₑ::MVector{D, Float64}
end

"""
    RotationParameters(dt, T, simulation, walkers,  H, rate, cage_time, escape_time)

General constructor for the `RotationParameters` structure.
"""
function RotationParameters(dt::Float64, T::Float64, simulation::String, Dᵣ::Float64, walkers::Int, H::Float64, rate::Float64, cage_time::Float64, escape_time::Float64)
    RotationParameters(dt, T, walkers, simulation, Dᵣ, H, rate, cage_time, escape_time)
end

"""
    RotationParameters(dt, T, walkers, order)

Constructor for a simple Brownian simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int; Dᵣ::Float64=1.0)
    tᵪ = MVector{walkers, Float64}(zeros(walkers))
    tₑ = MVector{walkers, Float64}([Inf for _ in 1:walkers])
    RotationParameters(dt, T, walkers,  "Brownian", Dᵣ, 0.0, 0.0, tᵪ, tₑ)
end

"""
    RotationParameters(dt, T, walkers,  H)

Constructor for a "Cage" simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int, H::Float64; Dᵣ::Float64=1.0)
    tᵪ = MVector{walkers, Float64}([Inf for _ in 1:walkers])
    tₑ = MVector{walkers, Float64}(zeros(walkers))
    RotationParameters(dt, T, walkers, "Cage", Dᵣ, H, 0.0, tᵪ, tₑ)
end

"""
    RotationParameters(dt, T, walkers,  H, rate)

Constructor for an "Escape" simulation.
"""
function RotationParameters(dt::Float64, T::Float64, walkers::Int, H::Float64, rate::Float64; Dᵣ::Float64=1.0)
    tᵪ = MVector{walkers, Float64}([Inf for _ in 1:walkers])
    tₑ = MVector{walkers, Float64}(zeros(walkers))
    RotationParameters(dt, T, walkers,  "Escape", Dᵣ, H, rate, tᵪ, tₑ)
end

# Utility Functions
include("io.jl")
include("math_utils.jl")
include("simulation.jl")
include("analysis.jl")

end # module