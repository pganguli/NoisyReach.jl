#!/usr/bin/env julia

begin
  import Pkg
  Pkg.activate("..")
  Pkg.instantiate()

  push!(LOAD_PATH, "$(@__DIR__)/../src")
  include("./pointcloud_data.jl")

  using NoisyReach
  using Random
  using Plots
  using Statistics

  Random.seed!(123) # Setting the seed
end

const sys = benchmarks[:F1]
const N = 500

const x0 = 1000.0
const u0 = 0.0

const T = 10
const hc = 1e-4
const Hc = floor(Int, T / hc)

experiment = ("casa", 10)

if length(ARGS) == 2
  experiment_name = ARGS[1]
  experiment_value = parse(Int, ARGS[2])
  experiment = (experiment_name, experiment_value)
end

hd = pointcloud_data[experiment[1]][experiment[2]]["h"]
Hd = floor(Int, T / hd)

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
xc[end], uc[end]

sys_d, K_d = synthesize(sys, hd)

z₀ = fill(x0, size(sys.A, 1))
xds, uds = unzip([simulate(sys_d, Hd, z₀, K_uncertain, K_d, pointcloud_data[experiment[1]][experiment[2]]["rel_err"]) for _ in 1:N])
xds[1][end], uds[1][end]

#rel_dev_trajs = [[(1 .- (xds[i][k] ./ xc[1+(k-1)*floor(Int,hd/hc)])) for k in 1:Hd] for i in 1:N]
dev_trajs = [[(xc[1+(k-1)*floor(Int,hd/hc)] .- xds[i][k]) for k in 1:Hd] for i in 1:N]
dev₁_trajs = [[dev[1] for dev in dev_traj] for dev_traj in dev_trajs]
uds₁ = [[udi[1] for udi in ud] for ud in uds]

mean_mean_dev₁ = mean([mean(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])
max_max_dev₁ = maximum([maximum(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])

mean_rmse_dev₁ = mean([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
std_rmse_dev₁ = std([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
max_rmse_dev₁ = maximum([std(dev₁_traj) for dev₁_traj in dev₁_trajs])

mean_mean_u₁ = mean([mean(abs.(ud₁)) for ud₁ in uds₁])
std_mean_u₁ = std([mean(abs.(ud₁)) for ud₁ in uds₁])
max_max_u₁ = maximum([maximum(abs.(ud₁)) for ud₁ in uds₁])

println("$(experiment[1]) @ $(experiment[2])")
println("mean_rmse_dev₁	std_rmse_dev₁	mean_mean_u₁	std_mean_u₁")
println("$mean_rmse_dev₁	$std_rmse_dev₁	$mean_mean_u₁	$std_mean_u₁")

convergences = [first_convergence([0.0, 0.0], xds[i], threshold=1e-2) for i in 1:N]
#median(convergences)

#plot_trajectories(xds[8], xc, hd, hc, T)

dev₁s = [dev₁ for dev₁ in dev₁_trajs[8]]
#plot(1:Hd, dev₁s, label="\$\\Delta x_1\$", lw=2)

#println("Trajectories that failed to converge:\n", [i for i in 1:N if isnothing(first_convergence([0.0, 0.0], xds[i], threshold=1e-2))])

#show(stdout, "text/plain", uds[5])

#maximum(xd -> maximum(xd -> abs(xd[1]), xd), xds)

#maximum(ud -> maximum(ud -> abs(ud[1]), ud), uds)

#maximum(xc -> abs(xc[1]), xc)

#maximum(u -> abs(u[1]), uc)
