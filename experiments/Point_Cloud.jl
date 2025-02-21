begin
  import Pkg
  Pkg.activate("..")
  Pkg.instantiate()

  push!(LOAD_PATH, "$(@__DIR__)/../src")
  include("./pointcloud_data.jl")

  using NoisyReach
  using Random
  using Plots

  Random.seed!(123) # Setting the seed
  pythonplot()
end

experiment = ("casa", 0)

const sys = benchmarks[:F1]
const N = 500

const x0 = 1000.0
const u0 = 0.0

const T = 10
const hc = 1e-4
const hd = pointcloud_data[experiment[1]][experiment[2]]["h"]
const Hc = floor(Int, T / hc)
const Hd = floor(Int, T / hd)

@debug """-----BEGIN PARAMS-----
sys=$sys
N=$N
x0=$x0
u0=$u0
T=$T
hc=$hc
hd=$hd
Hc=$Hc
Hd=$Hd
-----END PARAMS-----"""

sys_c, K_c = synthesize(sys, hc)
x₀ = fill(x0, size(sys.A, 1))
xc, uc = simulate(sys_c, Hc, x₀, K_certain, K_c)

sys_d, K_d = synthesize(sys, hd)
z₀ = fill(x0, size(sys.A, 1))
xds, uds = unzip([simulate(sys_d, Hd, z₀, K_uncertain, K_d, pointcloud_data[experiment[1]][experiment[2]]["rel_err"]) for _ in 1:N])

println("Trajectories that failed to converge:\n", [i for i in 1:N if isnothing(first_convergence([0.0, 0.0], xds[i], threshold=1e-2))])

plot(0:hc:T, hcat(xc...)', lab=["x₁(t)" "x₂(t)"], xlabel="Time [s]", linecolor=[:red :blue], alpha=0.6)
plot!(0:hd:T, hcat(xds[1]...)', lab=["x₁[k]" "x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, marker=:circle, markersize=2)
gui()
