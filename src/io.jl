function initialize_trajectory!(filename, params::RotationParameters)
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
    end

end

function get_time_trajectory(filename)
    T, dt = load_param(filename, "T"), load_param(filename, "dt")
    N = Int(T / dt)
    return collect(dt:dt:N*dt)
end

function save_timestep!(filename::String, M::Vector{<:AbstractMatrix}, time)
    
    v = [prod(euler_from_rotation(m)) for m in M]
    jldopen(filename, "a") do f
        f["TimeSteps/$(time)"] = v
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

function load_param(filename, key)
    out = nothing
    jldopen(filename) do f
        g = f["Params"]
        out = g[key]
    end
    return out
end