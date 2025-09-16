using MoonWalk
using Random, Distributions, LinearAlgebra, StaticArrays, LaTeXStrings, Measures, JLD2

walkers = 100
T = 1e1
dt= 0.01

params = RotationParameters(dt, T, walkers)
simulation(params; path="examples/brownian/dt$(dt)T$(T)")