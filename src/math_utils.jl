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

function sample_exponential(params::RotationParameters; shift::Float64=0.0)
    if params.rate == 0.0
        params.cage_time = shift + 0.0

        params.escape_time = Inf
    elseif params.rate == Inf
        params.cage_time =  Inf
        params.escape_time = shift + 0.0
    else
        params.cage_time = shift + rand(Exponential(params.rate))
        params.escape_time = params.cage_time + params.dt
    end
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
    return θ, e
end

function is_ϕ_equal_θ_plus_2kpi_e(ϕ::AbstractVector, θ::Real, e::AbstractVector; tol=1e-8)
    k = round(dot(ϕ-θ*e, e) / 2π)
    return k
end

function log_spaced_interval(start_val, stop_val, num_points)
    interval = exp10.(range(log10(start_val), log10(stop_val), length=num_points))
    interval = [round(Int, elt) for elt in interval]
    return unique(interval)
end

function estimate_diffusion_coefficient(t::Vector, msd::Vector)
    @assert length(t) == length(msd) "Time and MSD vectors must have the same length"
    fit_range = msd .> 1

    t_fit = t[fit_range]
    msd_fit = msd[fit_range]
    log_fit_range = log_spaced_interval(1, length(t_fit), 10)
    t_fit = t_fit[log_fit_range]
    msd_fit = msd_fit[log_fit_range]
    A = hcat(t_fit, ones(length(t_fit)))  # Design matrix for linear regression
    coeffs = A \ msd_fit  # Solve for coefficients

    return coeffs[1] / 6  # D = slope / 2
end

function find_best_dϕ(ϕ, θ_test, θ, e_test, e, k)
    k_test = nothing
    best_m = nothing
    min_norm = Inf
    best_dϕ = nothing
    dΩ = θ_test*e_test - θ*e
    dθ = θ_test - θ
    best_dϕ = nothing
    
    if abs(dθ) < π
        k_test = k
    elseif dθ > 0
        k_test = k - 1
    else
        k_test = k + 1

    end
    #=
    for m in -abs(k)-3:abs(k)+3
        candidate = dΩ +  2π * (m * e_test - k * e)

        nrm = norm(candidate)
        if nrm < min_norm
            min_norm = nrm
            best_dϕ = dΩ + 2π * (m * e_test - k * e)
            best_m = m
            
        end
    end
    =#
    #best_dϕ = dΩ 
    best_dϕ = dΩ + 2π * (k_test * e_test - k * e)
    norm_diff = norm(best_dϕ)
    if norm_diff > π && e != zeros(3)
        return error("e_test $e_test and e $e are too far apart (norm diff $norm_diff ), cannot find best dϕ")
    end 
    return best_dϕ, abs(dθ) > π
end
