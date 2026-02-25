using MoonWalk

walkers = 1000
T = 10.0
dt= 0.001

params = RotationParameters(dt, T, walkers)
simulation(params; path="brownian")
