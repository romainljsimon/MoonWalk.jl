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