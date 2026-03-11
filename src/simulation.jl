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


    clock = 1
    while clock <= params.T

        if (clock ∈ scheduler)
            save_timestep(trajectory_file, angle_definitions, clock)
        end


        if isa(params, BrownianParameters)
            # Every step, move of a random amount
            dΩ = params.amplitude * randn(rng, Float64, 3)
            dR = rotation_matrix_from_omega(dΩ)
            R = R * dR
            clock += 1
        end

        for definition in angle_definitions
            step!(definition, R, dR)
        end

    end

    return nothing
end
