using MoonWalk
using ProgressMeter

T = 10000.0
amplitude = 0.1
nb_walkers = 1000

prog = Progress(nb_walkers; desc="Simulating walkers...")
for i in 1:nb_walkers
    params = BrownianParameters(T, amplitude)
    simulation(params; path="brownian/$i")
    next!(prog)
end
