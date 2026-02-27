function cage(Rₜ, dΩ, H)
    dR = rotation_matrix_from_omega(dΩ)
    θₜ, _ = euler_from_rotation(Rₜ * dR)
    if abs(θₜ) > H
        dR = I(3)
    end
    return dR
end

function simulation(params::RotationParameters; path::String="./", rng=Xoshiro(), number_of_points_per_decade::Int=10)
    mkpath(path)
    trajectory_file = joinpath(path, "traj.csv")

    N = round(Int, params.T / params.dt)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]
    Rₜ = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]
    dR = [SMatrix{3,3,Float64}(I) for _ in 1:params.walkers]

    number_of_decades = floor(log10(N) + 1)
    logrange(1, N, Int(number_of_decades * number_of_points_per_decade))
    scheduler = unique([round(x) for x in logrange(1, N, Int(number_of_decades * number_of_points_per_decade))])

    angle_definitions = [ExactRotation(params.walkers), IntegralDefinition(params.walkers), UnboundedDefinition(params.walkers)]

    n_point = 100
    angle_distirbutions = [AngleDistribution(n_point) for _ in 1:length(angle_definitions)]

    initialize_trajectory!(trajectory_file, params, angle_definitions)
    save_timestep!(trajectory_file, R, dR, angle_definitions, 0, scheduler, params.dt)

    if params.simulation == "Escape"
        sample_exponential!(params; rng=rng)
    end

    i = 1
    prog = Progress(N; desc="Simulating walkers...")
    while i < N+1

        dΩ = sqrt(params.dt * params.Dᵣ) *randn(rng, Float64, (3, params.walkers))

        if params.simulation == "Brownian"
            for walker in 1:params.walkers
                dR[walker] = rotation_matrix_from_omega(dΩ[:, walker])
                R[walker] = R[walker] * dR[walker]
            end
        elseif params.simulation == "Cage"
            for walker in 1:params.walkers
                dR[walker] = cage(R[walker], dΩ[:, walker], params.H)
                R[walker] = R[walker] * dR[walker]
            end
        elseif params.simulation == "Escape"
            # Rt stores the rotation matrix since the last jump
            # This defines the cage: we refuse any Rt that would be greater than H
            # It's reset to identity at each jump
            for walker in 1:params.walkers
                if i * params.dt < params.tᵪ[walker]
                    dR[walker] = cage(Rₜ[walker], dΩ[:, walker], params.H)
                elseif i * params.dt < params.tₑ[walker]
                    dR[walker] = rotation_matrix_from_omega(dΩ[:, walker] ./ (5 * sqrt(params.dt)))
                else
                    Rₜ[walker] = I(3)
                    sample_exponential!(params; rng=rng, shift=i*params.dt, i=[walker])
                end
                R[walker] = R[walker] * dR[walker]
                Rₜ[walker] = Rₜ[walker] * dR[walker]
            end
        elseif params.simulation == "CTRW"
            for walker in 1:params.walkers
                # Jump
                if i * params.dt > params.tₑ[walker]
                    dR[walker] = rotation_matrix_from_omega(dΩ[:, walker])
                    R[walker] = R[walker] * dR[walker]
                    sample_power_law_jump!(params; rng=rng, shift=i*params.dt, i=[walker])
                end
            end
        else
            error("Unknown simulation type")
        end

        save_timestep!(trajectory_file, R, dR, angle_definitions, i, scheduler, params.dt)

        # Accumulate statistics about distributions at regular intervals
        if i != 0 && i % 100 == 0
            for (definition, distribution) in zip(angle_definitions, angle_distirbutions)
                add_omegas!(distribution, get_omegas(definition))
            end
        end

        i += 1
        next!(prog)
    end

    # Write distributions
    for (definition, distribution) in zip(angle_definitions, angle_distirbutions)
        file_path = joinpath(path, "angle_$(definition.name).csv")
        write_distribution(distribution, file_path)
    end

    return nothing
end
