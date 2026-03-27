using MoonWalk
using ProgressMeter

T = 100000.0
amplitude_small = 0.05
amplitude_large = 0.5
αs = [0.5, 0.8, 1.5]
cage_size = 0.2
nb_walkers = 100
τ = 1000.0

for α in αs
    println("α = $α")
    prog = Progress(nb_walkers; desc="Simulating walkers...")
    for i in 1:nb_walkers
        params = ParetoParameters(T, amplitude_small, amplitude_large, α, cage_size, τ)
        simulation(params; path="pareto/$α/$i")
        next!(prog)
    end
end
