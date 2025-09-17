abstract type TrajectoryMethod end

function analyze_trajectory(method::TrajectoryMethod, filename)
    T, dt = load_param(filename, "T"), load_param(filename, "dt")
    N = Int(T / dt)
    ϕ = Vector{Any}(undef, N)  # preallocate results

    for i in 1:N
        euler_vectors = load_timestep(filename, i)
        ϕ[i] = step!(method, euler_vectors)
    end
    return ϕ
end

struct ThetaMethod <: TrajectoryMethod
    n::Vector{Vector{Float64}}  # state carried across timesteps
end

function ThetaMethod(num_vectors::Int)
    # initialize n with unit vectors along x-axis
    return ThetaMethod([[1, 0, 0] for _ in 1:num_vectors])
end

function step!(method::ThetaMethod, euler_vectors)
    ϕᵢ = Vector{Float64}(undef, length(euler_vectors))
    for (i, ev) in enumerate(euler_vectors)
        θ = norm(ev)
        n = ev / θ
        pm = sign(dot(n, method.n[i]))
        ϕᵢ[i] = θ * pm
        method.n[i] = n * pm
    end
    return ϕᵢ
end

mutable struct IntegralMethod <: TrajectoryMethod
    Ωarray::Vector{Vector{Float64}}  # state carried across timesteps
    ϕarray::Vector{Vector{Float64}}
end

function IntegralMethod(num_vectors::Int)
    # initialize n with unit vectors along x-axis
    return IntegralMethod([[0, 0, 0] for _ in 1:num_vectors], [[0, 0, 0] for _ in 1:num_vectors])
end

function step!(method::IntegralMethod, Ωarray)
    ϕᵢ = copy(method.ϕarray)
    for (i, Ω) in enumerate(Ωarray)
        Rₜ = rotation_matrix_from_omega(Ω)
        R = rotation_matrix_from_omega(method.Ωarray[i])
        dR = transpose(R) * Rₜ
        dθ, n =  euler_from_rotation(dR)
        dϕ = dθ * n
        ϕᵢ[i] += dϕ 
    end
    method.ϕarray = ϕᵢ
    method.Ωarray = Ωarray
    return ϕᵢ
end