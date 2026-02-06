function initialize_trajectory!(filename, params::RotationParameters, scheduler)
    file_handle = jldopen(filename, "w")

    JLD2.Group(file_handle, "TimeSteps")

    m = JLD2.Group(file_handle, "Params")
    m["walkers"] = params.walkers
    m["H"] = params.H
    m["rate"] = params.rate
    m["tᵪ"] = params.tᵪ
    m["tₑ"] = params.tₑ
    m["dt"] = params.dt
    m["T"] = params.T
    m["scheduler"] = scheduler

    return file_handle
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

function save_timestep!(file_handle::JLD2.JLDFile{JLD2.MmapIO}, M::Vector{<:AbstractMatrix}, time, scheduler)
    if (time ∈ scheduler) || (time == 0)
        v = [prod(euler_from_rotation(m)) for m in M]
        file_handle["TimeSteps/$(time)"] = v
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


function read_lammpstrj_by_molecule(filename)
    positions = Dict{Int, Dict{Int, Dict{String, Any}}}() 
    current_timestep = nothing
    reading_atoms = false
    atom_data = Float64[]
    ncols = 9  # id type molecule x y z ix iy iz

    open(filename, "r") do io
        for line in eachline(io)
            if startswith(line, "ITEM: TIMESTEP")
                # lire le numéro du timestep
                current_timestep = parse(Int, readline(io))

            elseif startswith(line, "ITEM: ATOMS")
                # on commence à lire les lignes d'atomes
                reading_atoms = true
                atom_data = Float64[]

            elseif startswith(line, "ITEM:") && reading_atoms
                # fin de la section ATOMS : transformer atom_data en matrice
                natoms = length(atom_data) ÷ ncols
                data_matrix = reshape(atom_data, ncols, natoms)'  # N_atoms × 6
                xyz = data_matrix[:, 4:6]        # colonnes x y z
                types = round.(Int, data_matrix[:, 2])
                mol_ids = Int.(data_matrix[:, 3])
                ixiyiz = round.(Int, data_matrix[:, 7:9])

                # Regrouper par molécule
                mol_dict = Dict{Int, Dict{String, Any}}()
                for mol in unique(mol_ids)
                    mask = mol_ids .== mol
                    mol_dict[mol] = Dict(
                        "xyz" => xyz[mask, :],
                        "type" => types[mask],
                        "image" => ixiyiz[mask, :]
                    )
                end

                positions[current_timestep] = mol_dict
                reading_atoms = false

            elseif reading_atoms
                # lire une ligne d'atome et ajouter au vecteur
                append!(atom_data, parse.(Float64, split(line)))
            end
        end

        # gérer le dernier bloc ATOMS si le fichier se termine après
        if reading_atoms && current_timestep !== nothing
            natoms = length(atom_data) ÷ ncols
            data_matrix = reshape(atom_data, ncols, natoms)'
            xyz = data_matrix[:, 4:6]
            types = round.(Int, data_matrix[:, 2])
            mol_ids = Int.(data_matrix[:, 3])
            ixiyiz = round.(Int, data_matrix[:, 7:9])

            mol_dict = Dict{Int, Dict{String, Any}}()
            for mol in unique(mol_ids)
                mask = mol_ids .== mol
                mol_dict[mol] = Dict(
                    "xyz" => xyz[mask, :],
                    "type" => types[mask],
                    "image" => ixiyiz[mask, :]
                )
            end

            positions[current_timestep] = mol_dict
        end
    end

    return positions
end
