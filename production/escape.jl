using MoonWalk
using ProgressMeter

T = 100000.0
amplitude_small = 0.1
amplitude_large = 0.05
rate = 500.0
cage_size = 0.2
nb_walkers = 1000

prog = Progress(nb_walkers; desc="Simulating walkers...")
for i in 1:nb_walkers
    params = CageEscapeParameters(T, amplitude_small, amplitude_large, rate, cage_size)
    simulation(params; path="escape/$i")
    next!(prog)
end
