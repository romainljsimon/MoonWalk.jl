abstract type TrajectoryMethod end

function analyze_trajectory(method::TrajectoryMethod, filename; scheduler=1)
    T, dt = load_param(filename, "T"), load_param(filename, "dt")
    N = Int(T / dt)
    scheduler = set_scheduler(scheduler, N)
    t, ϕ = [], []
    for i in 1:N
        Ωarray = load_timestep(filename, i)
        result = step!(method, Ωarray)
        if i ∈ scheduler
            push!(ϕ, result)
            push!(t, i * dt)
        end
    end
    return t, ϕ
end

function create_log_scheduler(N)
    logspaced = 10 .^(range(0,stop=log10(N),length=50))
    scheduler = [Int(floor(x)) for x in logspaced]
    unique!(scheduler)
    return scheduler
end

function set_scheduler(x::Int, N)
    out = collect(1:x:N)
    return out
end

function set_scheduler(x::Vector{Int}, N)
    sort!(x)
    max_scheduler = maximum(x)
    shift = max_scheduler
    out = copy(x)
    while shift + max_scheduler < N
        out = vcat(out, shift .+ init_scheduler)
        shift += max_scheduler
    end
    return out
end


struct ThetaMethod <: TrajectoryMethod
    n::Vector{Vector{Float64}}  # state carried across timesteps
end

function ThetaMethod(num_vectors::Int)
    # initialize n with unit vectors along x-axis
    return ThetaMethod([[1, 0, 0] for _ in 1:num_vectors])
end

function step!(method::ThetaMethod, Ωarray)
    ϕᵢ = Vector{Float64}(undef, length(Ωarray))
    for (i, Ω) in enumerate(Ωarray)
        θ = norm(Ω)
        n = Ω / θ
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