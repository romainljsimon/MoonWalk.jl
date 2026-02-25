abstract type AngleDefinition end


struct ExactRotation <: AngleDefinition
    R::Vector{SMatrix{3,3,Float64}}
    name::String
end

function ExactRotation(n_walker::Int)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:n_walker]
    return ExactRotation(R, "ExactRotation")
end

function step!(method::ExactRotation, Rs::Vector{<:AbstractMatrix}, dRs::Vector{<:AbstractMatrix})
    for (i, R) in enumerate(Rs)
        method.R[i] = R
    end
end

function get_omegas(method::ExactRotation)
    return [prod(euler_from_rotation(m)) for m in method.R]
end

struct IntegralDefinition <: AngleDefinition
    R::Vector{SMatrix{3,3,Float64}}
    ϕs::Vector{Vector{Float64}}
    name::String
end

function IntegralDefinition(n_walker::Int)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:n_walker]
    ϕs = [[0, 0, 0] for _ in 1:n_walker]
    return IntegralDefinition(R, ϕs, "Integral")
end

function step!(method::IntegralDefinition, Rs::Vector{<:AbstractMatrix}, dRs::Vector{<:AbstractMatrix})
    for (i, R) in enumerate(Rs)
        dR = dRs[i]
        dθ, n =  euler_from_rotation(dR)
        dϕ = dθ * n
        method.ϕs[i] += dϕ
        method.R[i] = R
    end
end

function get_omegas(method::IntegralDefinition)
    return method.ϕs
end

struct UnboundedDefinition <: AngleDefinition
    R::Vector{SMatrix{3,3,Float64}}
    ϕs::Vector{Vector{Float64}}
    ns::Vector{Vector{Float64}}
    θs::Vector{Float64}
    name::String
end

function UnboundedDefinition(n_walker::Int)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:n_walker]
    ϕs = [[0, 0, 0] for _ in 1:n_walker]
    ns = [[1, 0, 0] for _ in 1:n_walker]
    θs = [0.0 for _ in 1:n_walker]
    return UnboundedDefinition(R, ϕs, ns, θs, "Unbounded")
end

function step!(method::UnboundedDefinition, Rs::Vector{<:AbstractMatrix}, dRs::Vector{<:AbstractMatrix})
    for (i, R) in enumerate(Rs)
        # This is the rotation matrix between the current position
        # and an origin defined by R[i]
        dR = transpose(R) * method.R[i]
        θ, n =  euler_from_rotation(dR)
        # If I go above the threshold
        if θ > 3
            # Cumulate
            method.ϕs[i] += θ * n
            # And reset the origin
            method.R[i] = R
            method.θs[i] = 0
            method.ns[i] = n
        else
            method.θs[i] = θ
            method.ns[i] = n
        end
    end
end

function get_omegas(method::UnboundedDefinition)
    ϕs_from_last_origin = method.θs .* method.ns
    return method.ϕs .+ ϕs_from_last_origin
end
