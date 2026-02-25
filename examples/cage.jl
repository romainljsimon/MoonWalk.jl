using MoonWalk

walkers = 1000
T = 10.0
dt= 0.001
H = 0.5

params = RotationParameters(dt, T, walkers, H)
simulation(params; path="cage")
