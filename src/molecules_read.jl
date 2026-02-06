function jld2_conv(filename,temp,p)
    traj = read_lammpstrj_by_molecule(filename)
    times = sort(collect(keys(traj)))[2:end]
    T = length(times)
    dt = 1
    N_mol = 1000
    path = "real/Temp$(temp)/dt$(dt)T$(T)_$(p)"
    mkpath(path)
    params = RotationParameters(Float64(dt), Float64(T), N_mol, Dᵣ=1.0)
    trajectory_file = joinpath(path, "traj.jld2") 
    file_handle = MoonWalk.initialize_trajectory!(trajectory_file, params, 1:T)
    R = [SMatrix{3,3,Float64}(I) for _ in 1:N_mol]
    for (i,j) in enumerate(times)
        for mol in 1:N_mol
            R[mol] = rot_from_xyz(traj,j,mol)
        end
        save_timestep!(file_handle, R, i,1:T)
    end
    return nothing
end




function tri(traj,timestep,molecule)
    box = 1.3572088e+01
    mol = traj[timestep][molecule]
    v=[(mol["xyz"][i,:],mol["type"][i]) for i in 1:3]
    sort!(v, by = x -> x[2])
    x₁,x₂,x₃ = v[1][1][1],v[2][1][1],v[3][1][1]
    y₁,y₂,y₃ = v[1][1][2],v[2][1][2],v[3][1][2]
    z₁,z₂,z₃ = v[1][1][3],v[2][1][3],v[3][1][3]
    if abs(x₁-x₂)>box/2
        x₂ += sign(x₁-x₂)*box
    end
    if abs(x₁-x₃)>box/2
        x₃ += sign(x₁-x₃)*box
    end   
    if abs(y₁-y₂)>box/2
        y₂ += sign(y₁-y₂)*box
    end
    if abs(y₁-y₃)>box/2
        y₃ += sign(y₁-y₃)*box
    end   
    if abs(z₁-z₂)>box/2
        z₂ += sign(z₁-z₂)*box
    end
    if abs(z₁-z₃)>box/2
        z₃ += sign(z₁-z₃)*box
    end   
    A = [x₁,y₁,z₁]
    B = [x₂,y₂,z₂]
    C = [x₃,y₃,z₃]
    return A,B,C
end


function tri_pos(traj,timestep,molecule)
    box = 1.3572088e+01
    mol = traj[timestep][molecule]
    v=[(mol["xyz"][i,:],mol["type"][i],mol["image"][i,:]) for i in 1:3]
    sort!(v, by = x -> x[2])
    for i in 1:3
        v[i][1] .+= box*v[i][3]
    end
    x₁,x₂,x₃ = v[1][1][1],v[2][1][1],v[3][1][1]
    y₁,y₂,y₃ = v[1][1][2],v[2][1][2],v[3][1][2]
    z₁,z₂,z₃ = v[1][1][3],v[2][1][3],v[3][1][3]
    A = [x₁,y₁,z₁]
    B = [x₂,y₂,z₂]
    C = [x₃,y₃,z₃]
    return A,B,C
end



function normal_vec(A,B,C)
    v1 = A .- C
    v2 = B .- C
    n = cross(v1,v2)/norm(cross(v1,v2))
    return n
end


function theta_from_mol(n1,n2)
    ps = dot(n1,n2)
    if abs(ps) > 1
        ps = sign(ps)
    end
    θ = acos(ps)
    n = cross(n1,n2)/norm(cross(n1,n2))
    return θ,n
end


function rot_from_xyz(traj,timestep,molecule)
    mol_1 = tri(traj,1,molecule)
    mol_t = tri(traj,timestep,molecule)
    θ,n = theta_from_mol(normal_vec(mol_1[1],mol_1[2],mol_1[3]),normal_vec(mol_t[1],mol_t[2],mol_t[3]))
    Ω = θ * n
    R = rotation_matrix_from_omega(Ω)
    return R
end

function theta_from_xyz(traj,timestep,molecule)
    mol_1 = tri(traj,1,molecule)
    mol_t = tri(traj,timestep,molecule)
    θ,n = theta_from_mol(normal_vec(mol_1[1],mol_1[2],mol_1[3]),normal_vec(mol_t[1],mol_t[2],mol_t[3]))
    Ω = θ * n
    R = rotation_matrix_from_omega(Ω)
    return Ω
end

