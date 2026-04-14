using MoonWalk
using Random

# Free diffusion
T = 1000.0
amplitude_small = 0.1

params = BrownianParameters(T, amplitude_small)
rng = Xoshiro(43)
simulation(params; path="trajectory/brownian/", rng=rng, save_positions=true)

# Cage
cage_size = 0.2
params = CageParameters(T, amplitude_small, cage_size)
rng = Xoshiro(42)
simulation(params; path="trajectory/cage/", rng=rng, save_positions=true)

# Cage escape
rate = 400.0
amplitude_large = 1
T = 3000.0

rng = Xoshiro(54) #   50 (with T=1300), 54 (with T=3000)
params = CageEscapeParameters(T, amplitude_small, amplitude_large, rate, cage_size)
simulation(params; path="trajectory/escape", rng=rng, save_positions=true)
