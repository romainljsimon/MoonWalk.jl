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


function save_timestep!(filename,Rs::Vector{<:AbstractMatrix},dRs::Vector{<:AbstractMatrix}, definitions::Vector{<:AngleDefinition}, time, scheduler)

    for definition in definitions
        step!(definition, Rs, dRs)
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
