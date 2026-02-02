abstract type AngleDefinition end


struct ExactRotation <: AngleDefinition
    R::Vector{SMatrix{3,3,Float64}}
    name::String
end

function ExactRotation(n_walker::Int)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:n_walker]
    return ExactRotation(R, "ExactRotation")
end

function step!(method::ExactRotation, M::Vector{<:AbstractMatrix})
    for (i, m) in enumerate(M)
        method.R[i] = m
    end
end

function get_omegas(method::ExactRotation)
    return [prod(euler_from_rotation(m)) for m in method.R]
end
