abstract type AngleDefinition end


mutable struct ExactRotation <: AngleDefinition
    R::SMatrix{3,3,Float64}
    name::String
end

function ExactRotation()
    R = SMatrix{3,3,Float64}(I)
    return ExactRotation(R, "ExactRotation")
end

function step!(method::ExactRotation, R::SMatrix, dR::SMatrix)
    method.R = R
end

function get_omega(method::ExactRotation)
    return prod(euler_from_rotation(method.R))
end

mutable struct IntegralDefinition <: AngleDefinition
    ϕ::SVector{3, Float64}
    name::String
end

function IntegralDefinition()
    ϕ = [0, 0, 0]
    return IntegralDefinition(ϕ, "Integral")
end

function step!(method::IntegralDefinition, R::SMatrix, dR::SMatrix)
    dθ, n =  euler_from_rotation(dR)
    dϕ = dθ * n
    method.ϕ += dϕ
end

function get_omega(method::IntegralDefinition)
    return method.ϕ
end

mutable struct UnboundedDefinition <: AngleDefinition
    R::SMatrix{3,3,Float64}
    ϕ::SVector{3, Float64}
    n::SVector{3, Float64}
    θ::Float64
    name::String
end

function UnboundedDefinition()
    R = SMatrix{3,3,Float64}(I)
    ϕ = [0, 0, 0]
    n = [1, 0, 0]
    θ = 0.0
    return UnboundedDefinition(R, ϕ, n, θ, "Unbounded")
end

function step!(method::UnboundedDefinition, R::SMatrix, dR::SMatrix)
    # This is the rotation matrix between the current position
    # and an origin defined by the method's R
    dR = transpose(R) * method.R
    θ, n =  euler_from_rotation(dR)
    # If I go above the threshold
    if θ > 3
        # Cumulate
        method.ϕ += θ * n
        # And reset the origin
        method.R = R
        method.θ = 0
        method.n = n
    else
        method.θ = θ
        method.n = n
    end
end

function get_omega(method::UnboundedDefinition)
    ϕ_from_last_origin = method.θ * method.n
    return method.ϕ + ϕ_from_last_origin
end
