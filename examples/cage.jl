using MoonWalk

walkers = 1000
T = 10.0
dt= 0.001

Hs = [0.5, 0.8, 1.0]

for H in Hs
    pritnln("H $H")
    params = RotationParameters(dt, T, walkers, H)
    simulation(params; path="examples/cage/t$(dt)T$(T)H$(H)")
end
