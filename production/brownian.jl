using MoonWalk
using Random

T = 100000.0
amplitude = 0.1
nb_walkers = 10000

i = parse(Int, ARGS[1])
println("Processing simulation $i out of $nb_walkers")

params = BrownianParameters(T, amplitude)
rng = Xoshiro(i)
simulation(params; path="brownian/$i", rng=rng)
