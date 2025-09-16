function initialize_trajectory!(filename, params::RotationParameters)
    jldopen(filename, "w") do f
        JLD2.Group(f, "TimeSteps")
        f["walkers"] = params.walkers
        f["H"] = params.H
        f["rate"] = params.rate
        f["tᵪ"] = params.tᵪ
        f["tₑ"] = params.tₑ
        f["dt"] = params.dt
        f["T"] = params.T
    end

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