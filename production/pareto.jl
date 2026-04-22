using MoonWalk
using Random

amplitude_small = 0.2
amplitude_large = 0.5
cage_size = 0.2
nb_walkers = 200000
τ = 10000.0

αs = [0.7, 1.2]

total_number_of_simulations = nb_walkers * length(αs)
i = parse(Int, ARGS[1])
walker_id = (i - 1) % nb_walkers + 1
α = αs[Int(floor((i - 1) / nb_walkers)) + 1]

if isfile("pareto/$α/$walker_id/traj.csv")
    exit()
end

if walker_id <= 2000
   T = 100000000.0
   number_of_points_in_output=100
elseif walker_id > 20000
   T = 1000000.0
   number_of_points_in_output=78
else
   T = 10000000.0
   number_of_points_in_output=89
end

println("Processing simulation $i out of $total_number_of_simulations")
println("Walker n°$walker_id, α = $α")

params = ParetoParameters(T, amplitude_small, amplitude_large, α, cage_size, τ)
rng = Xoshiro(i)
simulation(params; path="pareto/$α/$walker_id", rng=rng, number_of_points_in_output=number_of_points_in_output)
