using MoonWalk
using Random

T = 100000000.0
amplitude_small = 0.05
amplitude_large = 0.5
cage_size = 0.2
nb_walkers = 10000
τ = 100000.0

αs = [0.7, 1.2]

total_number_of_simulations = nb_walkers * length(αs)
i = parse(Int, ARGS[1])
walker_id = (i - 1) % nb_walkers + 1
α = αs[Int(floor((i - 1) / nb_walkers)) + 1]

println("Processing simulation $i out of $total_number_of_simulations")
println("Walker n°$walker_id, α = $α")

params = ParetoParameters(T, amplitude_small, amplitude_large, α, cage_size, τ)
rng = Xoshiro(i)
simulation(params; path="pareto/$α/$walker_id", rng=rng)
