
struct AngleDistribution
    dx::Float64
    bin_centers::Vector{Float64}
    bin_counters::Vector{Int64}
end

function AngleDistribution(n_bin::Int)
    dx = 2π / n_bin
    bin_centers = [-π + dx * (i - 0.5) for i in 1:n_bin]
    bin_counters = [0 for _ in 1:n_bin]
    return AngleDistribution(dx, bin_centers, bin_counters)
end


function add_omegas!(distrib::AngleDistribution, omegas::Vector{SVector{3, Float64}})
    for omega in omegas
        for component in omega
            idx = Int(ceil(mod(component + π, 2π) / distrib.dx))
            distrib.bin_counters[idx] += 1
        end
    end
end

function write_distribution(distrib::AngleDistribution, filepath::String)
    open(filepath, "w") do file
        write(file, "theta,p\n")

        integral = distrib.dx * sum(distrib.bin_counters)
        for (theta, counter) in zip(distrib.bin_centers, distrib.bin_counters)
            value = counter / integral
            write(file, "$(theta),$(value)\n")
        end
    end
end
