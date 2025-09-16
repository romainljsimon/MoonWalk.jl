using MoonWalk

walkers = 100
T = 1e1
dt= 0.01

params = RotationParameters(dt, T, walkers)
simulation(params; path="examples/brownian/dt$(dt)T$(T)")