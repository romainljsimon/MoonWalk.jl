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
