
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

function unwrap_home(angles::AbstractVector{T}; range = π) where T
    unwrapped = copy(angles)
    for i in 2:length(angles)
        diff = unwrapped[i] - unwrapped[i-1]
        if diff > range
            unwrapped[i:end] .+= -diff
        elseif diff < -range
            unwrapped[i:end] .-= diff
        end
    end
    return unwrapped
end

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

function bch_new(ϕ, wiener, order)
    corr  = zero(ϕ)
    skew_ϕ = skew_symmetric(ϕ)
    skew_wiener = skew_symmetric(wiener)
    if order >= 1
        pw = skew_ϕ * wiener
        corr += 1.0 / 2.0 * pw
    end
    if order >= 2
        ppw = skew_ϕ * skew_ϕ * wiener
        wpw = skew_wiener * skew_ϕ * wiener
        corr += 1.0/12.0*(ppw - wpw)
    end
    if order >= 3
        wppw= skew_wiener * skew_ϕ * skew_ϕ * wiener
        corr += -1.0/24.0 * wppw
    end
    if order >= 4
        pxpxpxpxw = - dot(ϕ, ϕ) * pxpxw
        wxwxwxwxp = dot(wiener, wiener) * wxpxw
        corr += -1.0/720.0 * (pxpxpxpxw + wxwxwxwxp)
        wxpxpxpxw = - dot(ϕ, ϕ) * wxpxw
        pxwxwxwxp = dot(wiener, wiener) * pxpxw
        corr += 1.0/360.0 * (wxpxpxpxw + pxwxwxwxp)
        wxpxwxpxw =  - dot(ϕ, wiener) * wxpxw
        pxwxpxwxp = dot(wiener, ϕ) * pxpxw
        corr += 1.0 /120 * (wxpxwxpxw + pxwxpxwxp)
    end
    if any(corr .> 1000) || any(corr .< -1000)
        return NaN, NaN, NaN
    else
        return corr
    end
    return corr
end

function is_ϕ_equal_θ_plus_2kpi_e(ϕ::AbstractVector, θ::Real, e::AbstractVector; tol=1e-8)
    k = round(dot(ϕ-θ*e, e) / 2π)
    return k
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
    # Retourner le vecteur propre correspondant
    e = real(eigenvectors[:, idx])
    sin_θ = - tr(skew_symmetric(e)*R)/2
    θ = atan(sin_θ, cos_θ)
    #θ = acos(clamp((tr(R) - 1)/ 2, -1.0, 1.0) )
    if isapprox(θ, 0.0; atol=1e-8)
        return 0.0, SVector(0.0, 0.0, 0.0)
    end
    # Calculer l'axe de rotation normalisé
    #e = SVector(R[3,2] - R[2,3], R[1,3] - R[3,1], R[2,1] - R[1,2]) / (2 * sin(θ))
    return θ, e
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

function cage(ϕ, R, R_thresh, wiener, H, order)
    dR = rotation_matrix_from_omega(wiener)
    R_test = R * dR
    angle_test, _ = euler_from_rotation(R_thresh * dR)
    if abs(angle_test) < H
        corr = bch_new(ϕ, wiener, order)
        corrnan = isnan(corr[1]) || isnan(corr[2]) || isnan(corr[3])
        if corrnan
            dϕ = NaN, NaN, NaN
        else
            dϕ = wiener + corr
        end
    else
        dϕ = 0
        dR = I(3)
    end
    return dϕ, dR
end

function cage(ϕ, R, R_thresh, wiener, H, ::InfinityInteger)
    dR = rotation_matrix_from_omega(wiener)
   
    θ_tresh, e_thresh = euler_from_rotation(R_thresh * dR)
    if abs(θ_tresh) < H
        R_test = R * dR
        θ_test, e_test = euler_from_rotation(R_test)
        θ, e = euler_from_rotation(R)
        k = is_ϕ_equal_θ_plus_2kpi_e(ϕ, θ, e)
        θ_test, e_test = euler_from_rotation(R_test)
        if dot(e_test, e) < 0
            e_test = -e_test
            θ_test = -θ_test
        end
        dϕ = find_best_dϕ(ϕ, θ_test, θ, e_test, e, k)
    else
        dϕ = 0
        dR = I(3)
    end
    return dϕ, dR
end


function brownian(ϕ, R, wiener, order)
    dR = rotation_matrix_from_omega(wiener)
    R_test = R * dR
    corr = bch_new(ϕ, wiener, order)
    corrnan = isnan(corr[1]) || isnan(corr[2]) || isnan(corr[3])
    if corrnan
        dϕ = NaN, NaN, NaN
    else
        
        dϕ = wiener + corr
        R = copy(R_test)
    end
    return dϕ, dR
end

function find_best_dϕ(ϕ, θ_test, θ, e_test, e, k)
    k_test = nothing
    best_m = nothing
    min_norm = Inf
    best_dϕ = nothing
    dΩ = θ_test*e_test - θ*e
    dθ = θ_test - θ
    best_dϕ = nothing
    

    for m in -abs(k)-3:abs(k)+3
        candidate = ϕ + dΩ +  2π * (m * e_test - k * e)

        nrm = abs(norm(candidate) - norm(ϕ))
        if nrm < min_norm
            min_norm = nrm
            best_dϕ = dΩ + 2π * (m * e_test - k * e)
            best_m = m
            
        end
    end
    return best_dϕ
end

function brownian(ϕ, R, wiener, ::InfinityInteger)
    
    dR = rotation_matrix_from_omega(wiener)
    R_test = R * dR
    θ, e = euler_from_rotation(R)
    k = is_ϕ_equal_θ_plus_2kpi_e(ϕ, θ, e)
    θ_test, e_test = euler_from_rotation(R_test)
    if dot(e_test, e) < 0
        e_test = -e_test
        θ_test = -θ_test
    end
    dϕ = find_best_dϕ(ϕ, θ_test, θ, e_test, e, k)
    return dϕ, dR
end

function solve_sde(params::RotationParameters, seed=42)
    #Random.seed!(seed)
    ϕ_0 = zeros(3)
    N = round(Int, params.T / params.dt)  # Number of time steps
    t = range(0, stop=params.T , length=N)
    msd = zeros(N)
    walkers = params.walkers
    dϕ_matrix = zeros(walkers, N)  # Matrice pour stocker dϕ pour chaque walker
    x_matrix, y_matrix, z_matrix = zeros(walkers, N), zeros(walkers, N), zeros(walkers, N)
    frobenius = zeros(N)
    
    for walker in 1:walkers
        println("Walker ", walker)
        ϕ_matrix = zeros(N, 3)
        ϕ = copy(ϕ_0)
        R = I(3)
        R_thresh = I(3)
        dϕ = 0
        dR = I(3)
        dW = sqrt(params.dt) * randn(3, N)  # Wiener increments
        if params.simulation == "Escape"
            sample_exponential(params)
        end
        i = 2
        while i < N+1
            if i*params.dt < params.cage_time
                dϕ, dR = cage(ϕ, R, R_thresh, dW[:, i], params.H, params.order)
            elseif i*params.dt < params.escape_time
                if params.simulation == "Escape"
                    dϕ, dR = brownian(ϕ, R, dW[:, i] ./ (5 * sqrt(params.dt)), params.order)
                else
                    dϕ, dR = brownian(ϕ, R, dW[:, i], params.order)
                end
                
            else
                R_thresh = I(3)
                sample_exponential(params; shift=i*params.dt)
                continue
            end
            ϕ .+= dϕ
            R = R*dR
            R_thresh = R_thresh*dR
            frobenius_i =  sqrt(sum((rotation_matrix_from_omega(ϕ) - R).^2))
            frobenius[i] += frobenius_i / walkers
            ϕ_matrix[i, :] = ϕ
            i += 1
        end
        x_matrix[walker, :], y_matrix[walker, :], z_matrix[walker, :] = ϕ_matrix[:, 1], ϕ_matrix[:, 2], ϕ_matrix[:, 3]
        dϕ_matrix[walker, :] .= ϕ_matrix[:, 1] #(x[1:end].^2 + y[1:end].^2 + z[1:end].^2).^(0.5)
    end
    for (q, l)  in enumerate(eachcol(x_matrix))
        elt_x = mean(x_matrix[:, q].^2)
        elt_y = mean(y_matrix[:, q].^2)
        elt_z = mean(z_matrix[:, q].^2)
        msd[q] =  sum([elt_x, elt_y, elt_z])
    end
    return t, msd, (x_matrix, y_matrix, z_matrix), frobenius
end