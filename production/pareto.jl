using MoonWalk
using ProgressMeter

T = 10000000.0
amplitude_small = 0.05
amplitude_large = 0.5
cage_size = 0.1
nb_walkers = 1000
τ = 10000.0

αs = [1.5]

prog = Progress(nb_walkers * length(αs); desc="Simulating walkers...")
for α in αs
    Threads.@threads for i in 1:nb_walkers
        params = ParetoParameters(T, amplitude_small, amplitude_large, α, cage_size, τ)
        simulation(params; path="pareto/$α/$i")
        next!(prog)
    end
end
