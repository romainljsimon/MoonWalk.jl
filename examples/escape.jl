using MoonWalk

walkers = 100
H = 0.8
T = 1e1
dt= 0.01
rate = 1.0

params = RotationParameters(dt, T, walkers, H, rate)
simulation(params; path="examples/escape/dt$(dt)T$(T)H$(H)rate$(rate)")