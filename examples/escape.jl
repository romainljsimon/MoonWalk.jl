using MoonWalk

walkers = 100
H = 0.8
T = 100.0
dt= 0.001
rate = 2.0

params = RotationParameters(dt, T, walkers, H, rate)
simulation(params; path="escape")
