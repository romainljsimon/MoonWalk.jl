using MoonWalk
using ProgressMeter

T = 10000.0
amplitude = 0.1
cage_size = 1.0
nb_walkers = 10000

prog = Progress(nb_walkers; desc="Simulating walkers...")
for i in 1:nb_walkers
    params = CageParameters(T, amplitude, cage_size)
    simulation(params; path="cage/$i")
    next!(prog)
end
