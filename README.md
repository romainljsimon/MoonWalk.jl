# MoonWalk.jl

MoonWalk.jl is a Julia package for simulating stochastic differential equations (SDEs) on rotation groups. It provides tools for modeling Brownian motion, cage dynamics, and escape dynamics.

## Installation

To install the package, clone the repository then precompile it:

```bash
git clone https://github.com/romainljsimon/MoonWalk.jl
```

## Publication

This package was used to produce the results of ADD-CITATION-AND-LINK.

To reproduce the results, you need to run each julia script in the [production](./production) subfolder, which will run the diffusion simulations.
All the simulations are seeded, so running these scripts will produce the exact same data as in the paper.

Once this is done, you can run each python script in the [analysis](./analysis) subfolder to post-process the trajectories, compute MSDs, and create the corresponding plots.

Each subfolder contains an additional README file with more details.
