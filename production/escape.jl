using MoonWalk
using Random

T = 10000000.0
amplitude_small = 0.2
amplitude_large = 0.5
cage_size = 0.2
nb_walkers = 10000

rates = [3000000.0, 1000000.0, 300000.0, 100000.0, 30000.0, 10000.0, 3000.0, 1000.0, 300.0, 100.0]

total_number_of_simulations = nb_walkers * length(rates)
i = parse(Int, ARGS[1])
walker_id = (i - 1) % nb_walkers + 1
rate = rates[Int(floor((i - 1) / nb_walkers)) + 1]

println("Processing simulation $i out of $total_number_of_simulations")
println("Walker n°$walker_id, rate = $rate")

params = CageEscapeParameters(T, amplitude_small, amplitude_large, rate, cage_size)
rng = Xoshiro(i)
simulation(params; path="escape/$rate/$walker_id", rng=rng)
