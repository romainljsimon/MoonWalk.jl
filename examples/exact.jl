using MoonWalk
using Random, Distributions, Plots, LinearAlgebra, StaticArrays, LaTeXStrings, Measures
# Example usage
walkers = 10
H = 0.2
T = 1e2
dts= [0.01 ]
rate = 1.0
order = InfinityInteger()   # Order of BCH
min_dt = minimum(dts)
colors = [:blue, :red, :green, :purple, :orange, :black]

selected_times = [min_dt, T-1.0]  # Indices of the time steps (adjust based on your N)
traj = plot(layout=(1,1), size=(600, 300), framestyle=:box, grid=:false, tick_direction=:in)
fig = plot(layout=(1,2), framestyle=:box, grid=:false, legend=:topleft, 
tick_direction=:in, size=(600, 300), xlim = (min_dt,T), xscale=:log10, yscale=:log10, 
foreground_color_legend = nothing, left_margin = [5mm 5mm], bottom_margin = [3mm 0mm])
for (i, dt) in enumerate(dts)
    println("dt $dt")
    params = RotationParameters(dt, T, walkers, order)
    t, msd, tuple_matrix, frobenius = solve_sde(params)
    plot!(fig[1], t[2:end], msd[2:end], label=nothing)
    plot!(fig[2], t[2:end], frobenius[2:end], legend=:bottomright, label=nothing)
    plot!(traj, t, tuple_matrix[1]'[:, 1], xlabel="t", label="x")
    plot!(traj, t, tuple_matrix[2]'[:, 1], xlabel="t", label="y")
    plot!(traj, t, tuple_matrix[3]'[:, 1], xlabel="t", label="z")
end


#plot!()
plot!(fig[1], [min_dt, T], 3*[min_dt, T]/5,  ylim=(min_dt, 1e2), color=:black, linestyle=:dash, label=nothing)
plot!(fig[2], [min_dt, T], [2.40, 2.40], ylim=(1e-16, 1e1), color=:black, linestyle=:dash, label=nothing)
annotate!(fig[2], 1.7e-3, 4e-3 , Plots.text(L"d_F(e^{\tilde{\phi}(t)}, R(t))", rotation=90))
annotate!(fig[1], 2.3e-3, 3.3*1e-1, Plots.text(L"\langle \phi^2(t) \rangle", rotation=90))
annotate!(fig[2], 3e-2, 1, L"d_{rand}")
annotate!(fig[1], 0.35, 6.2e-3, L"t")
annotate!(fig[2], 0.35, 2e-17, L"t")
annotate!(fig[1], 0.05, 0.02, L"\propto t")

#x = range(-π,π)
#y = cos.(x)./4       
#plot!(hist_fig, x, y, xlabel="dϕ", ylabel="PDF", subplot=0)#, yscale=:log10)
#plot!(fig[1], [dt, T], 3*[dt, T]./rate, color=:black, linestyle=:dash, label="Brownian")
#t, msd, dϕ_matrix, frobenius = solve_sde(dt, T, 0, rate=rate, walkers=100)
#plot!(fig[1], [100*dt, T], ([100*dt, T]/20).^(1.75), color=:black, linestyle=:dash)
#display(hist_fig)
display(traj)
display(fig)
#savefig(fig, "./escapecage.pdf")