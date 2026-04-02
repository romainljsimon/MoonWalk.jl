using MoonWalk
using ProgressMeter

T = 1000000.0
amplitude_small = 0.05
amplitude_large = 0.1
cage_size = 0.1
nb_walkers = 100

rates = [300000.0, 100000.0, 30000.0, 10000.0, 3000.0, 1000.0, 300.0, 100.0, 30.0, 10.0]

prog = Progress(nb_walkers * length(rates); desc="Simulating walkers...")
Threads.@threads for rate in rates
    for i in 1:nb_walkers
        params = CageEscapeParameters(T, amplitude_small, amplitude_large, rate, cage_size)
        simulation(params; path="escape/$rate/$i")
        next!(prog)
    end
end
