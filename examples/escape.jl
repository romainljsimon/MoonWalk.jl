using MoonWalk

walkers = 1000
H = 0.8
T = 10.0
dt= 0.001

rates = [0.1, 1.0, 10.0]
for rate in rates
    println("Rate $rate")
    params = RotationParameters(dt, T, walkers, H, rate)
    simulation(params; path="examples/escape/dt$(dt)T$(T)H$(H)rate$(rate)")
end
