# Convert a 3D vector into a skew-symmetric matrix
function skew_symmetric(v::AbstractVector{T}) where T
    x, y, z = v
    return @SMatrix [
        0  -z   y
        z   0  -x
       -y   x   0
    ]
end


function rotation_matrix_from_omega(omega::AbstractVector{T}) where T
    θ = norm(omega)
    if θ < 1e-8
        return I(3)
    end
    axis = omega / θ
    ω_tilde = skew_symmetric(axis)

    # Rodrigues' formula
    R = I(3) + sin(θ) * ω_tilde + (1 - cos(θ)) * (ω_tilde * ω_tilde)
    return R
end


function omega_from_skew(R::AbstractMatrix)
    return SVector(R[3,2], R[1,3], R[2,1])
end


function euler_from_rotation(R::AbstractMatrix)
    cos_θ = clamp((tr(R) - 1)/ 2, -1.0, 1.0)
    eigen_decomp = eigen(R)
    eigenvalues = eigen_decomp.values
    eigenvectors = eigen_decomp.vectors
    idx = findfirst(x -> isapprox(x, 1.0; atol=1e-8), eigenvalues)
    @assert idx !== nothing "Aucune valeur propre égale à 1 n'a été trouvée"
    e = real(eigenvectors[:, idx])
    sin_θ = - tr(skew_symmetric(e)*R)/2
    θ = atan(sin_θ, cos_θ)
    if isapprox(θ, 0.0; atol=1e-8)
        return 0.0, SVector(0.0, 0.0, 0.0)
    end
    if θ > 0
        return θ, e
    else
        return -θ, -e
    end
end


function log_spaced_interval(start_val, stop_val, num_points)
    interval = exp10.(range(log10(start_val), log10(stop_val), length=num_points))
    return interval
end


function euler_angles_from_rotation(R::AbstractMatrix)
    theta = -asin(R[3,1])
    psi = atan(R[3,2]/cos(theta),R[3,3]/cos(theta))
    phi = atan(R[2,1]/cos(theta),R[1,1]/cos(theta))
    return psi,theta,phi
end


function smooth_pareto(rng::Xoshiro, params::ParetoParameters)

    flat_contribution = params.α / (1 + params.α)

    u = rand(rng)
    if u < flat_contribution
        return u / flat_contribution * params.τ
    else
        return rand(rng, Pareto(params.α, params.τ))
    end
end
