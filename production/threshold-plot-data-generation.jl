using MoonWalk
using LinearAlgebra, Random, StaticArrays

seed = parse(Int, ARGS[1])

max_time = 1000.0
rng = Xoshiro(seed)
amplitude = 0.5

unbounded_definition = UnboundedDefinition()
exact_definition = ExactRotation()

time = 0.0
R = SMatrix{3,3,Float64}(I)

file = open("angles.csv", "w")
write(file, "time,theta,theta_rebuilt,theta_local\n")

while time < max_time
    global time += MoonWalk.small_jump_distribution(rng)

    dΩ = amplitude * (rand(rng, Float64, 3) .- 0.5)
    dR = MoonWalk.rotation_matrix_from_omega(dΩ)
    global R = R * dR

    MoonWalk.step!(unbounded_definition, R, dR)
    MoonWalk.step!(exact_definition, R, dR)

    exact_θ = norm(MoonWalk.get_omega(exact_definition))
    unbounded_θ = norm(MoonWalk.get_omega(unbounded_definition))
    truncated_θ = unbounded_definition.θ
    write(file, "$time,$exact_θ,$unbounded_θ,$truncated_θ\n")
end
