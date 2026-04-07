function cage(Rₜ, dΩ, H)
    dR = rotation_matrix_from_omega(dΩ)
    θₜ, _ = euler_from_rotation(Rₜ * dR)
    if abs(θₜ) > H
        dR = I(3)
    end
    return dR
end

function simulation(params::RotationParameters; path::String="./", rng::Xoshiro=Xoshiro(), number_of_points_in_output::Int=100, save_positions::Bool=false)

    mkpath(path)
    trajectory_file = joinpath(path, "traj.csv")

    scheduler = log_spaced_interval(0.1, params.T, number_of_points_in_output)

    angle_definitions = [ExactRotation(), IntegralDefinition(), UnboundedDefinition()]


    initialize_trajectory!(trajectory_file, params, angle_definitions)

    if save_positions
        position_files = joinpath(path, "positions.csv")
        open(position_files, "w") do file
            write(file, "time,x,y,z\n")
        end
    end

    R = SMatrix{3,3,Float64}(I)
    dR = SMatrix{3,3,Float64}(I)

    # Need to initialize for cage escape
    if isa(params, CageEscapeParameters) || isa(params, ParetoParameters)
        time_of_next_small_jump = small_jump_distribution(rng)
        if isa(params, CageEscapeParameters)
            time_of_next_big_jump = rand(rng, Exponential(params.rate))
        else
            time_of_next_big_jump = smooth_pareto(rng, params)
        end
        Rₜ = SMatrix{3,3,Float64}(I)
    end


    clock = 0
    next_index_of_scheduler_for_printing = 1
    while clock <= params.T

        # Pure diffusive
        if isa(params, BrownianParameters)

            time_of_next_jump = clock + small_jump_distribution(rng)

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

            time_of_next_jump = clock + small_jump_distribution(rng)

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

        # Cage + big jumps (either according to exponential or Pareto distribution)
        elseif isa(params, CageEscapeParameters) || isa(params, ParetoParameters)

            time_of_next_jump = min(time_of_next_small_jump, time_of_next_big_jump)

            while time_of_next_jump > scheduler[next_index_of_scheduler_for_printing]
                save_timestep(trajectory_file, angle_definitions, scheduler[next_index_of_scheduler_for_printing])

                if next_index_of_scheduler_for_printing == length(scheduler)
                    break
                end

                next_index_of_scheduler_for_printing += 1
            end

            # Stay in cage
            if time_of_next_small_jump < time_of_next_big_jump
                dΩ = params.amplitude_small * (rand(rng, Float64, 3) .- 0.5)
                dR = rotation_matrix_from_omega(dΩ)

                # How much we moved since the origin
                θ, _ = euler_from_rotation(Rₜ * dR)
                if abs(θ) > params.cage_size
                    dR = SMatrix{3,3,Float64}(I)
                end

                clock = time_of_next_small_jump
                time_of_next_small_jump += small_jump_distribution(rng)

                Rₜ = Rₜ * dR
            # Larger jump
            else
                dΩ = params.amplitude_large * (rand(rng, Float64, 3) .- 0.5)
                dR = rotation_matrix_from_omega(dΩ)

                # Reset origin for the cage
                Rₜ = SMatrix{3,3,Float64}(I)

                clock = time_of_next_big_jump
                time_of_next_small_jump = clock + small_jump_distribution(rng)

                if isa(params, CageEscapeParameters)
                    time_of_next_big_jump = clock + rand(rng, Exponential(params.rate))
                else
                    time_of_next_big_jump = clock + smooth_pareto(rng, params)
                end

            end

            R = R * dR
        end


        for definition in angle_definitions
            step!(definition, R, dR)
        end

        if save_positions
            xyz = R * [1, 0, 0]
            open(position_files, "a") do file
                write(file, "$clock,$(xyz[1]),$(xyz[2]),$(xyz[3])\n")
            end
        end

    end

    return nothing
end
