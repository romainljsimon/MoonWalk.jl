using MoonWalk
using ProgressMeter

T = 100000.0
amplitude = 0.5
αs = [0.5, 0.8, 1.5]
nb_walkers = 1000

for α in αs
    println("α = $α")
    prog = Progress(nb_walkers; desc="Simulating walkers...")
    for i in 1:nb_walkers
        params = ParetoParameters(T, amplitude, α)
        simulation(params; path="pareto/$α/$i")
        next!(prog)
    end
end
