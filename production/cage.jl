using MoonWalk
using Random

T = 100000.0
amplitude = 0.2
cage_size = 0.2
nb_walkers = 10000

i = parse(Int, ARGS[1])
println("Processing simulation $i out of $nb_walkers")

params = CageParameters(T, amplitude, cage_size)
rng = Xoshiro(i)
simulation(params; path="cage/$i", rng=rng)
