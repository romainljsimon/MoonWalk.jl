function initialize_trajectory!(filename, params::RotationParameters, definitions::Vector{<:AngleDefinition})
    open(filename, "w") do file
        write(file, "time")
        for definition in definitions
            write(file, ",$(definition.name)_x")
            write(file, ",$(definition.name)_y")
            write(file, ",$(definition.name)_z")
        end
        write(file, "\n")
    end
end


function save_timestep(filename,definitions::Vector{<:AngleDefinition}, time)

    open(filename, "a") do file
        write(file, "$(time)")
        for definition in definitions
            vs = get_omega(definition)
            write(file, ",$(vs[1])")
            write(file, ",$(vs[2])")
            write(file, ",$(vs[3])")
        end
        write(file, "\n")
    end
end
