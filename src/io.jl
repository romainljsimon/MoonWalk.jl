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

function get_length_simulation(filename)
    T, dt = load_param(filename, "T"), load_param(filename, "dt")
    N = Int(T / dt)
    return N
end

function get_time_trajectory(filename)
    N = get_length_simulation(filename)
    return collect(dt:dt:N*dt)
end


function save_timestep!(filename,M::Vector{<:AbstractMatrix}, definitions::Vector{<:AngleDefinition}, time, scheduler)

    for definition in definitions
        step!(definition, M)
    end

    if (time ∈ scheduler) || (time == 0)
        open(filename, "a") do file
            write(file, "time")
            for definition in definitions
                v = get_omegas(definition)
                write(file, ",$(v[1])")
            end
            write(file, "\n")
        end
    end
end

function load_timestep(filename, time)
    out = nothing
    jldopen(filename) do f
        g = f["TimeSteps"]
        out = g["$(time)"]
    end
    return out
end

function load_trajectory_walker(filename, walker)
    N = get_length_simulation(filename)
    Ωwalker = []
    for i in i:N
        Ωarray = load_timestep(filename, i)
        push!(Ωwalker, Ωarray[walker])
    end
    return Ωwalker
end

function load_param(filename, key)
    out = nothing
    jldopen(filename) do f
        g = f["Params"]
        out = g[key]
    end
    return out
end

function load_timesteps(filename)
    return load_param(filename, "scheduler")
end
