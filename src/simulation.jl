function cage(Rₜ, dΩ, H)
    dR = rotation_matrix_from_omega(dΩ)
    θₜ, _ = euler_from_rotation(Rₜ * dR)
    if abs(θₜ) > H
        dR = I(3)
    end
    return dR
end

function simulation(params::RotationParameters; path::String="./", rng=Xoshiro())
    mkpath(path)
    trajectory_file = joinpath(path, "traj.jld2")
    
    N = round(Int, params.T / params.dt)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]
    Rₜ = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]
    dR = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]
    i = 1
    initialize_trajectory!(trajectory_file, params)
    save_timestep!(trajectory_file, R, 0)
    if params.simulation == "Escape"
        sample_exponential!(params; rng=rng)
    end
    prog = Progress(params.walkers; desc="Simulating walkers...")
    while i < N+1

        dΩ = sqrt(params.dt) * randn(rng, Float64, (3, params.walkers))
        for walker in 1:params.walkers
            if i * params.dt < params.tᵪ[walker]
                dR[walker] = cage(Rₜ[walker], dΩ[:, walker], params.H)
            elseif i * params.dt < params.tₑ[walker] 
                if params.simulation == "Escape"
                    dR[walker] = rotation_matrix_from_omega(dΩ[:, walker] ./ (5 * sqrt(params.dt)))
                else
                    dR[walker] = rotation_matrix_from_omega(dΩ[:, walker])
                end
            else
                Rₜ[walker] = I(3)
                sample_exponential!(params; rng=rng, shift=i*params.dt, i=[walker])
                continue
            end
            R[walker] = R[walker] * dR[walker]
            Rₜ[walker] = Rₜ[walker] * dR[walker]
            
        end
        save_timestep!(trajectory_file, R, i)
        i += 1
        next!(prog)
    end
    return nothing
end