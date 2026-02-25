function initialize_trajectory!(filename, params::RotationParameters, definitions::Vector{<:AngleDefinition})
    open(filename, "w") do file
        write(file, "# walkers=$(params.walkers),H=$(params.H),rate=$(params.rate),T=$(params.T)\n")

        write(file, "time")
        for definition in definitions
            write(file, ",$(definition.name)")
            write(file, "\n")
        end
    end
end


function load_omegas(filename::String)
    file_handle = jldopen(filename)

    params = file_handle["Params"]

    timesteps = params["scheduler"]
    n_walker = params["walkers"]
    omegas = [[zeros(SVector{3}) for _ in 1:n_walker] for _ in 1:length(timesteps)]

    for (i, time) in enumerate(timesteps)
        v = file_handle["TimeSteps/ExactRotation/$(time)"]
        for j in 1:n_walker
            omegas[i][j] = v[j]
        end
    end

    close(file_handle)

    return omegas, timesteps

end


function save_timestep!(filename,M::Vector{<:AbstractMatrix}, definitions::Vector{<:AngleDefinition}, time, scheduler)

    for definition in definitions
        step!(definition, M)
    end

    if (time ∈ scheduler) || (time == 0)
        open(filename, "a") do file
            write(file, "$time")
            for definition in definitions
                vs = get_omegas(definition)
                average_theta = mean([norm(v) for v in vs])
                write(file, ",$average_theta")
            end
            write(file, "\n")
        end
    end
end
