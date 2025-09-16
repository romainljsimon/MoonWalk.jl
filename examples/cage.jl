using MoonWalk

walkers = 100
H = 0.8
T = 1e1
dt= 0.01


params = RotationParameters(dt, T, walkers, H)
path = "examples/cage/t$(dt)T$(T)H$(H)"
simulation(params; path=path)

print(load_timestep(joinpath(path, "traj.jld2"), 3))