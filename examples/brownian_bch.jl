using MoonWalk
using Random, Distributions, Plots, LinearAlgebra, StaticArrays, LaTeXStrings, Measures
# Example usage
walkers = 1000
T = 1e2
dt= 0.01
orders = [2, 3]   # Order of BCH
colors = [:blue, :red, :green, :purple, :orange, :black]

selected_times = [dt, T-1.0]  # Indices of the time steps (adjust based on your N)
hist_fig = plot(layout=(length(orders), 1), size=(300, 300*length(orders)), framestyle=:box, grid=:false, tick_direction=:in)
traj = plot(layout=(1,1), size=(800, 300), framestyle=:box, grid=:false, tick_direction=:in)
fig = plot(layout=(1,2), framestyle=:box, grid=:false, legend=:topleft, 
tick_direction=:in, size=(600, 300), xlim = (dt,T), xscale=:log10, yscale=:log10, 
foreground_color_legend = nothing, left_margin = [5mm 5mm], bottom_margin = [3mm 0mm])

for (i, order) in enumerate(orders)
    println("Order $order")
    params = RotationParameters(dt, T, walkers, order)
    t, msd, dϕ_matrix, frobenius = solve_sde(params)
    plot!(fig[1], t[2:end], msd[2:end], label="order $order")
    plot!(fig[2], t[2:end], frobenius[2:end], legend=:bottomright, label="order $order")
end

plot!(fig[1], [dt, 0.1], 3*[dt, 0.1],  ylim=(dt, 1e2), color=:black, linestyle=:dash, label=nothing)
plot!(fig[2], [dt, T], [2.40, 2.40], ylim=(1e-6, 1e1),color=:black, linestyle=:dash, label=nothing)
annotate!(fig[2], 1.7e-3, 4e-3 , Plots.text(L"d_F(e^{\tilde{\phi}(t)}, R(t))", rotation=90))
annotate!(fig[1], 2.3e-3, 3.3*1e-1, Plots.text(L"\langle \phi^2(t) \rangle", rotation=90))
annotate!(fig[2], 3e-2, 1, L"d_{rand}")
annotate!(fig[1], 0.35, 6.2e-3, L"t")
annotate!(fig[2], 0.35, 2e-7, L"t")
annotate!(fig[1], 0.05, 0.02, L"\propto t")
display(fig)
savefig(fig, "examples/brownian_bch.pdf")