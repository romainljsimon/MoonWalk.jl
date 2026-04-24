# Running the simulations used for the publication

This directory contains six scripts.

Four of them (brownian.jl, cage.jl, escape.jl, pareto.jl) simulate random walkers on sphere, in different setups:
 - brownian.jl: the walkers can freely diffuse on the sphere
 - cage.jl: the walkers are confined to a cage on the sphere
 - escape.jl: the walkers are confined to a cage until they make a large jump (after a time drawn for an exponential distribution) and get trapped in a new cage. They alternate between staying trapped in a cage, and escaping it with a large jump
 - pareto.jl: similar to escape.jl, but the escape times are drawn from a modified Pareto distribution


Each of these scripts is used to simulate thousands of independent walkers. They must be run with an integer as an argument, which is the ID of the walker to process.

For example, pareto.jl is setup to run 400,000 walkers (200,000 for each of the two values of $\alpha$ to be used), so the script should be run 400,000 times, each time with a different argument, which must take its values between 1 and 400,000.


The last two scripts (create-trajectories.jl and threshold-plot-data-generation.jl) run very short trajectories, which are used to create visual representation of the our method.
