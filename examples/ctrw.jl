using MoonWalk

walkers = 100
T = 1000.0
dt= 0.01
α = 0.8
Dᵣ = 1.0

params = RotationParametersCTRW(dt, T, walkers, α; Dᵣ=Dᵣ)
simulation(params; path="ctrw")
