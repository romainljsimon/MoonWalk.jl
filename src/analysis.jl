abstract type TrajectoryMethod end

function analyze_trajectory(method::TrajectoryMethod, filename; scheduler=1)
    T, dt = load_param(filename, "T"), load_param(filename, "dt")
    N = Int(T / dt)
    scheduler = set_scheduler(scheduler, N)
    timesteps = load_timesteps(filename)
    t, ϕ, dfrob = [], [], []
    for i in timesteps
        Ωarray = load_timestep(filename, i)
        result = step!(method, Ωarray)
        if i ∈ scheduler
            push!(ϕ, result)
            push!(t, i * dt)
            push!(dfrob, mean_frobenius_norm(method.ϕarray, Ωarray))
        end
    end
    return t, ϕ, dfrob
end

function mean_frobenius_norm(ϕarray, Ωarray)
    out = 0.0
    for  (ϕ, Ω) in zip(ϕarray, Ωarray)
        out += norm(rotation_matrix_from_omega(ϕ) - rotation_matrix_from_omega(Ω))
    end
    return out
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
    ϕarray::Vector{Vector{Float64}}
end

function ThetaMethod(num_vectors::Int)
    # initialize n with unit vectors along x-axis
    return ThetaMethod([[1, 0, 0] for _ in 1:num_vectors], [[0.0, 0.0, 0.0] for _ in 1:num_vectors])
end

function step!(method::ThetaMethod, Ωarray)
    θᵢ = Vector{Float64}(undef, length(Ωarray))
    for (i, Ω) in enumerate(Ωarray)
        θ = norm(Ω)
        n = Ω / θ
        pm = sign(dot(n, method.n[i]))
        θᵢ[i] = θ * pm
        method.n[i] = n * pm
        method.ϕarray[i] = θ*n
    end
    return θᵢ, n
end

abstract type BCHOrder end
struct order0 <: BCHOrder end
struct order1 <: BCHOrder end
struct order2 <: BCHOrder end
struct order3 <: BCHOrder end

mutable struct IntegralMethod{O<:BCHOrder} <: TrajectoryMethod
    Ωarray::Vector{Vector{Float64}}  # state carried across timesteps
    ϕarray::Vector{Vector{Float64}}
    order::O
end



function IntegralMethod(num_vectors::Int; order::BCHOrder=order0())
    # initialize n with unit vectors along x-axis
    return IntegralMethod([[0., 0., 0.] for _ in 1:num_vectors], [[0., 0., 0.] for _ in 1:num_vectors], order)
end

function step!(method::IntegralMethod, Ωarray)
    ϕᵢ = copy(method.ϕarray)
    for (i, Ω) in enumerate(Ωarray)
        Rₜ = rotation_matrix_from_omega(Ω)
        R = rotation_matrix_from_omega(method.Ωarray[i])
        dR = transpose(R) * Rₜ
        dθ, n =  euler_from_rotation(dR)
        dϕ = dθ * n
        ϕᵢ[i] += bch(ϕᵢ[i], dϕ, method.order)
    end
    method.ϕarray = ϕᵢ
    method.Ωarray = Ωarray
    return ϕᵢ
end


mutable struct UnboundedThetaMethod <: TrajectoryMethod
    n::Vector{Vector{Float64}}
    θarray::Vector{Float64}
    ψarray::Vector{Float64}
    saut::Vector{Int}
    ϕarray::Vector{Vector{Float64}}
end

function UnboundedThetaMethod(num_vectors::Int)
    # initialize n with unit vectors along x-axis
    return ThetaMethod2([[1, 0, 0] for _ in 1:num_vectors], [0 for _ in 1:num_vectors],[0 for _ in 1:num_vectors],[0 for _ in 1:num_vectors],[[0., 0., 0.] for _ in 1:num_vectors])
end


function step!(method::UnboundedThetaMethod, Ωarray)
    θᵢ = copy(method.θarray)
    ψᵢ = Vector{Float64}(undef, length(Ωarray))
    nᵢ = copy(method.n)
    sᵢ = copy(method.saut)
    for (j, Ω) in enumerate(Ωarray)
        θ = norm(Ω)
        n = Ω / θ
        pm = sign(dot(n, nᵢ[j]))
        if abs(θ*pm - θᵢ[j])>3
            sᵢ[j] += sign(θᵢ[j])
        end
        nᵢ[j] = n * pm
        θᵢ[j] = θ * pm
        ψᵢ[j] = 2*pi*sᵢ[j] + θᵢ[j]
    end
    method.θarray = θᵢ
    method.ψarray = ψᵢ
    method.n = nᵢ
    method.saut = sᵢ
    method.ϕarray = ψᵢ .* nᵢ 
    return ψᵢ,nᵢ
end

mutable struct EulerMethod <: TrajectoryMethod
    Ωarray::Vector{Vector{Float64}}
    Rarray::Vector{Matrix{Float64}}
    relance::Vector{Vector{Float64}}
    ϕarray::Vector{Vector{Float64}}
end

function EulerMethod(num_vectors::Int)
    return EulerMethod([[0, 0, 0] for _ in 1:num_vectors],[[1 0 0 ; 0 1 0 ; 0 0 1] for _ in 1:num_vectors],[[0.0 , 0.0 , 0.0] for _ in 1:num_vectors],[[0.0 , 0.0 , 0.0] for _ in 1:num_vectors])
end 

function step!(method::EulerMethod, Ωarray)
    ψᵢ = Vector{Float64}(undef, length(Ωarray))
    θᵢ = Vector{Float64}(undef, length(Ωarray))
    ϕᵢ = Vector{Float64}(undef, length(Ωarray))
    for (j, Ω) in enumerate(Ωarray)
        R = rotation_matrix_from_omega(method.Ωarray[j])
        Rₜ = rotation_matrix_from_omega(Ω)
        dR = transpose(R) * Rₜ
        method.Rarray[j] = method.Rarray[j] * dR
        C = euler_angles_from_rotation(method.Rarray[j])
        ψᵢ[j] = C[1] + method.relance[j][1]
        θᵢ[j] = C[2] + method.relance[j][2]
        ϕᵢ[j] = C[3] + method.relance[j][3]
        if abs(C[1]) > 1 || abs(C[2]) > 1 || abs(C[3]) > 1
            method.relance[j] .+= [C[1],C[2],C[3]]
            method.Rarray[j] = I(3) 
        end 
        method.ϕarray[j] = rotation_matrix_from_omega(method.Rarray[j])
    end
    method.Ωarray = Ωarray
    return ψᵢ ,θᵢ ,ϕᵢ
end
function bch(ϕᵢ, dϕ, ::order0)
    return dϕ
end

function bch(ϕᵢ, dϕ, ::order1)
    dϕₒ = bch(ϕᵢ, dϕ, order0())
    return dϕₒ + 0.5 * cross(ϕᵢ, dϕ)
end

function bch(ϕᵢ, dϕ, ::order2)
    dϕₒ = bch(ϕᵢ, dϕ, order1())
    cross1 = cross(ϕᵢ, dϕ)
    dϕₒ += 1. / 12. * cross(ϕᵢ, cross1)
    dϕₒ -= 1. / 12. * cross(dϕ, cross1)
    return dϕₒ
end

function bch(ϕᵢ, dϕ, ::order3)
    dϕₒ = bch(ϕᵢ, dϕ, order2())
    cross1 = cross(ϕᵢ, dϕ)
    cross2 = cross(ϕᵢ, cross1)
    return dϕₒ - 1. / 24. * cross(dϕ, cross2)
end



mutable struct ThreshThetaMethod{O<:BCHOrder} <: TrajectoryMethod
    Ωarray::Vector{Vector{Float64}}  # state carried across timesteps
    ϕarray::Vector{Vector{Float64}}
    ϕthresharray::Vector{Vector{Float64}}
    thresh::Float64
    order::O
end



function ThreshThetaMethod(num_vectors::Int; thresh::Float64=1.0, order::BCHOrder=order0())
    # initialize n with unit vectors along x-axis
    return ThreshThetaMethod([[0., 0., 0.] for _ in 1:num_vectors], [[0., 0., 0.] for _ in 1:num_vectors], [[0., 0., 0.] for _ in 1:num_vectors], thresh, order)
end

function step!(method::ThreshThetaMethod, Ωarray)
    ϕᵢ = copy(method.ϕthresharray)
    for (i, Ω) in enumerate(Ωarray)
        Rₜ = rotation_matrix_from_omega(Ω)
        R = rotation_matrix_from_omega(method.Ωarray[i])
        dR = transpose(R) * Rₜ
        dθ, n =  euler_from_rotation(dR)
        dϕ = dθ * n
        ϕᵢ[i] += bch(ϕᵢ[i], dϕ, method.order)
        if dθ > method.thresh
            method.Ωarray[i] = Ω
            method.ϕthresharray[i] += bch(method.ϕthresharray[i], dϕ, method.order)
        end
    end
    method.ϕarray = ϕᵢ
    return ϕᵢ
end
