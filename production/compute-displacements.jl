using MoonWalk
using DataStructures
using LinearAlgebra

filename = "examples/brownian/dt0.001T10.0/traj.jld2"

dt= 0.001
omegas, timesteps = load_omegas(filename)

ϕ_sq = zeros(Float64, length(timesteps))
nb_walker = length(omegas[1])
for (i, omega_walker) in enumerate(omegas)
    ϕ_sq[i] = sum(norm.(omega_walker).^2) / nb_walker
end

open("data-brownian.txt", "w") do io
    write(io, "#t phi\n")
    for (t, phi) in zip(timesteps, ϕ_sq)
        write(io, "$(t*dt) $phi\n")
    end
end
