function initialize_trajectory!(filename, params::RotationParameters, scheduler)
    jldopen(filename, "w") do f
        JLD2.Group(f, "TimeSteps")
        m = JLD2.Group(f, "Params")
        m["walkers"] = params.walkers
        m["H"] = params.H
        m["rate"] = params.rate
        m["tᵪ"] = params.tᵪ
        m["tₑ"] = params.tₑ
        m["dt"] = params.dt
        m["T"] = params.T
        m["scheduler"] = scheduler
    end

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

function save_timestep!(filename::String, M::Vector{<:AbstractMatrix}, time, scheduler)
    if (time ∈ scheduler) || (time == 0)
        v = [prod(euler_from_rotation(m)) for m in M]
        jldopen(filename, "a") do f
            f["TimeSteps/$(time)"] = v
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