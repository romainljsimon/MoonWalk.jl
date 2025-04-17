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
    if any(corr .> 1000) || any(corr .< -1000)
        return NaN, NaN, NaN
    else
        return corr
    end
    return corr
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
    ϕ_0 = zeros(3)
    N = round(Int, params.T / params.dt)
    t = range(0, stop=params.T , length=N)
    msd = zeros(N)
    walkers = params.walkers
    dϕ_matrix = zeros(walkers, N)
    x_matrix, y_matrix, z_matrix = zeros(walkers, N), zeros(walkers, N), zeros(walkers, N)
    frobenius = zeros(N)
    
    for walker in 1:params.walkers
        println("Walker ", walker)
        ϕ_matrix = zeros(N, 3)
        ϕ = copy(ϕ_0)
        R = I(3)
        R_thresh = I(3)
        dϕ = 0
        dR = I(3)
        dW = sqrt(params.dt) * randn(3, N)
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
        dϕ_matrix[walker, :] .= ϕ_matrix[:, 1]
    end
    for (q, l)  in enumerate(eachcol(x_matrix))
        elt_x = mean(x_matrix[:, q].^2)
        elt_y = mean(y_matrix[:, q].^2)
        elt_z = mean(z_matrix[:, q].^2)
        msd[q] =  sum([elt_x, elt_y, elt_z])
    end
    return t, msd, (x_matrix, y_matrix, z_matrix), frobenius
end