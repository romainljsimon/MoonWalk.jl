using MoonWalk
using ProgressMeter

T = 10000.0
amplitude = 0.1
cage_size = 1.0
nb_walkers = 10000

i = parse(Int, ARGS[1])
println("Processing simulation $i out of $nb_walkers")

params = CageParameters(T, amplitude, cage_size)
simulation(params; path="cage/$i")
