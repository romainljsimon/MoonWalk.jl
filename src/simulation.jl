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

    number_of_decades = floor(log10(params.T) + 1)
    scheduler = unique([round(x) for x in logrange(1, params.T, Int(number_of_decades * number_of_points_per_decade))])

    angle_definitions = [ExactRotation(), IntegralDefinition(), UnboundedDefinition()]


    initialize_trajectory!(trajectory_file, params, angle_definitions)


    R = SMatrix{3,3,Float64}(I)
    Rₜ = SMatrix{3,3,Float64}(I)
    dR = SMatrix{3,3,Float64}(I)


    clock = 0
    next_index_of_scheduler_for_printing = 1
    while clock <= params.T

        # Pure diffusive
        if isa(params, BrownianParameters)

            time_of_next_jump = clock + rand(rng)

            while time_of_next_jump > scheduler[next_index_of_scheduler_for_printing]
                save_timestep(trajectory_file, angle_definitions, scheduler[next_index_of_scheduler_for_printing])

                if next_index_of_scheduler_for_printing == length(scheduler)
                    break
                end

                next_index_of_scheduler_for_printing += 1
            end


            dΩ = params.amplitude * (rand(rng, Float64, 3) .- 0.5)
            dR = rotation_matrix_from_omega(dΩ)
            R = R * dR

            clock = time_of_next_jump

        # Trapped in a cage
        elseif isa(params, CageParameters)

            time_of_next_jump = clock + rand(rng)

            while time_of_next_jump > scheduler[next_index_of_scheduler_for_printing]
                save_timestep(trajectory_file, angle_definitions, scheduler[next_index_of_scheduler_for_printing])

                if next_index_of_scheduler_for_printing == length(scheduler)
                    break
                end

                next_index_of_scheduler_for_printing += 1
            end


            dΩ = params.amplitude * (rand(rng, Float64, 3) .- 0.5)
            dR = rotation_matrix_from_omega(dΩ)

            # How much we moved since the origin
            θ, _ = euler_from_rotation(R * dR)
            if abs(θ) > params.cage_size
                dR = SMatrix{3,3,Float64}(I)
            end

            R = R * dR

            clock = time_of_next_jump

        end

        for definition in angle_definitions
            step!(definition, R, dR)
        end

    end

    return nothing
end
