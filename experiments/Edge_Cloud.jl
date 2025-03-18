### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "Edge Cloud Controllers"
#> description = "Designing controllers for split edge-cloud networks"
#> 
#>     [[frontmatter.author]]
#>     name = "Prateek Ganguli"
#>     url = "https://pganguli.github.io"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6d3ba880-d5f7-4cdf-83db-fd321031f867
# ╠═╡ show_logs = false
begin
  import Pkg
  Pkg.activate("..")
  Pkg.instantiate()

  push!(LOAD_PATH, "$(@__DIR__)/../src")

  using NoisyReach
  using Random
  using Plots
  using Statistics

  Random.seed!(123) # Setting the seed

  using PlutoUI
end

# ╔═╡ 78de028c-1eb1-4120-a3c5-455754ffa9fb
md"""
# Import Packages
"""

# ╔═╡ 08bed08d-7ce8-4f07-b055-c1e5b9e33ade
html"""
<style>
  img {
    max-width: 100% !important;
    max-height: 100% !important;
  }
</style>
"""

# ╔═╡ cb10d46e-378a-4726-90ea-e1f4b4be796f
md"""
# Configure Experiment
"""

# ╔═╡ 3cd33b19-87de-4280-8a5f-73f071c08c4b
md"""
## Common Parameters
"""

# ╔═╡ ca5f1534-4201-4dab-a917-9d0d15b0f13d
begin
  sys = benchmarks[:F1]
  const N = 500
  const x0 = 1000.0
  const u0 = 0.0
  const hc = 1e-4
  nothing
end

# ╔═╡ f73455e6-e275-4b02-93fe-cb311a60fdc4
md"""
## Tunable Parameters
"""

# ╔═╡ 565778df-82b2-448d-a6bd-626d8fe7fc31
T = @bind T Slider(0.1:0.1:10, default=2.5, show_value=true)

# ╔═╡ 3635e9ac-8429-49bf-93dd-96eec97b807f
hd = @bind hd Slider(0:0.01:T, default=0.02, show_value=true)

# ╔═╡ 9e2f99ea-b711-492f-915d-4981325856e7
α = @bind α Slider(0:0.1:1, default=0.5, show_value=true)

# ╔═╡ ac159d66-9878-4aca-aa88-0dd576481010
Dcₑ = @bind Dcₑ Slider(0:0.001:hd, default=0.01, show_value=true)

# ╔═╡ 24abfc73-8ff4-4a2b-b377-78c134fea835
Dcₔ = @bind Dcₔ Slider(0:0.001:5, default=0.1, show_value=true)

# ╔═╡ f421fe0b-b121-4c74-9ea2-b45c328759b6
σₔ = @bind σₔ Slider(0:0.05:10, default=0.3, show_value=true)

# ╔═╡ 8f06880f-0046-4569-a790-19b03988f36a
σₑ = @bind σₑ Slider(σₔ:0.05:10, default=0.7, show_value=true)

# ╔═╡ 20d653a3-9489-4671-9432-ddad7e83bc8a
μ = @bind μ Slider(-5:0.05:5, default=0, show_value=true)

# ╔═╡ 0d3bfeca-b728-466a-a8c9-4b4187072d55
md"""
## Computed Parameters
"""

# ╔═╡ 2d2f7f38-fe61-44d4-b3a4-0bd584da7c1a
begin
  @show const n = floor(Int, Dcₔ / hd)
  @show const Hc = floor(Int, T / hc)
  @show const Hd = floor(Int, T / hd)
  nothing
end

# ╔═╡ e9b3e724-0703-4aec-a8be-1af5da0a9ebe
md"""
# Calculate Trajectory
"""

# ╔═╡ cb49ecee-ff53-4fa2-98be-e752885de5e2
md"""
## Ideal Trajectory
"""

# ╔═╡ 51add47f-1c67-4ac3-a3a2-d5d1f870c1a7
begin
  sys_c, K_c = synthesize(sys, hc)

  x₀ = fill(x0, size(sys.A, 1))
  xc, uc = simulate(sys_c, Hc, x₀, K_certain, K_c)
  nothing
end

# ╔═╡ 672d6769-0077-4c3a-9991-7f70b72edac4
md"""
## Only Edge Input
"""

# ╔═╡ c2cca18a-1434-464f-b774-06f217d00271
let
begin
  z₀ = [fill(x0, size(sys.A, 1)); u0]
  sys_d, K_d = synthesize(sys, hd, Dcₑ)

  xds, uds = unzip([simulate(sys_d, Hd, z₀, K_uncertain, K_d, σₑ, μ) for _ in 1:N])
  xds = [[xd_i[1:2] for xd_i in xd] for xd in xds]
  uds = [[ud_i[1] for ud_i in ud] for ud in uds]
end

begin
  #rel_dev_trajs = [[(1 .- (xds[i][k] ./ xc[1+(k-1)*floor(Int,hd/hc)])) for k in 1:Hd] for i in 1:N]
  dev_trajs = [[(xc[1+(k-1)*floor(Int, hd / hc)] .- xds[i][k]) for k in 1:Hd] for i in 1:N]
  dev₁_trajs = [[dev[1] for dev in dev_traj] for dev_traj in dev_trajs]
  uds₁ = [[udi[1] for udi in ud] for ud in uds]

  max_max_x₁ = maximum(xd -> maximum(xd -> abs(xd[1]), xd), xds)

  mean_mean_dev₁ = mean([mean(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])
  max_max_dev₁ = maximum([maximum(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])

  mean_rmse_dev₁ = mean([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  std_rmse_dev₁ = std([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  max_rmse_dev₁ = maximum([std(dev₁_traj) for dev₁_traj in dev₁_trajs])

  mean_mean_u₁ = mean([mean(abs.(ud₁)) for ud₁ in uds₁])
  std_mean_u₁ = std([mean(abs.(ud₁)) for ud₁ in uds₁])
  max_max_u₁ = maximum([maximum(abs.(ud₁)) for ud₁ in uds₁])

  convergences = [first_convergence([0.0, 0.0], xds[i], threshold=1e-2) for i in 1:N]
  mean_convergence_time = mean(filter(!isnothing, convergences)) * hd

  max_xc₁ = maximum(xc -> abs(xc[1]), xc)
  max_uc₁ = maximum(u -> abs(u[1]), uc)
end

begin
  println("mean_rmse_dev₁	std_rmse_dev₁	mean_mean_u₁	std_mean_u₁")
  println("$mean_rmse_dev₁	$std_rmse_dev₁	$mean_mean_u₁	$std_mean_u₁")
  println("")
  println("max_max_x₁	max_xc₁	max_max_u₁	max_uc₁")
  println("$max_max_x₁	$max_xc₁	$max_max_u₁	$max_uc₁")
  println("")
  println("mean_convergence_time")
  println("$mean_convergence_time")
 
  plot_trajectories(xds, xc, hd, hc, T), plot(1:Hd, dev₁_trajs[1], lab="\$\\Delta x_1\$", lw=2)
end
end

# ╔═╡ 86014118-222d-4cbe-b4ce-b90c6e6d9078
md"""
## Only Cloud Input
"""

# ╔═╡ 68b22a80-d58b-47a5-ad16-83591913839b
let
begin
  z₀ = [fill(x0, size(sys.A, 1) * (n + 1)); u0]
  sys_d, K_d = synthesize(sys, hd, Dcₔ)

  xds, uds = unzip([simulate(sys_d, Hd, z₀, K_uncertain, K_d, σₔ, μ) for _ in 1:N])
  xds = [[xd_i[1:2] for xd_i in xd] for xd in xds]
  uds = [[ud_i[1] for ud_i in ud] for ud in uds]
end

begin
  #rel_dev_trajs = [[(1 .- (xds[i][k] ./ xc[1+(k-1)*floor(Int,hd/hc)])) for k in 1:Hd] for i in 1:N]
  dev_trajs = [[(xc[1+(k-1)*floor(Int, hd / hc)] .- xds[i][k]) for k in 1:Hd] for i in 1:N]
  dev₁_trajs = [[dev[1] for dev in dev_traj] for dev_traj in dev_trajs]
  uds₁ = [[udi[1] for udi in ud] for ud in uds]

  max_max_x₁ = maximum(xd -> maximum(xd -> abs(xd[1]), xd), xds)

  mean_mean_dev₁ = mean([mean(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])
  max_max_dev₁ = maximum([maximum(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])

  mean_rmse_dev₁ = mean([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  std_rmse_dev₁ = std([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  max_rmse_dev₁ = maximum([std(dev₁_traj) for dev₁_traj in dev₁_trajs])

  mean_mean_u₁ = mean([mean(abs.(ud₁)) for ud₁ in uds₁])
  std_mean_u₁ = std([mean(abs.(ud₁)) for ud₁ in uds₁])
  max_max_u₁ = maximum([maximum(abs.(ud₁)) for ud₁ in uds₁])

  convergences = [first_convergence([0.0, 0.0], xds[i], threshold=1e-2) for i in 1:N]
  mean_convergence_time = mean(filter(!isnothing, convergences)) * hd
  
  max_xc₁ = maximum(xc -> abs(xc[1]), xc)
  max_uc₁ = maximum(u -> abs(u[1]), uc)
end

begin
  println("mean_rmse_dev₁	std_rmse_dev₁	mean_mean_u₁	std_mean_u₁")
  println("$mean_rmse_dev₁	$std_rmse_dev₁	$mean_mean_u₁	$std_mean_u₁")
  println("")
  println("max_max_x₁	max_xc₁	max_max_u₁	max_uc₁")
  println("$max_max_x₁	$max_xc₁	$max_max_u₁	$max_uc₁")
  println("")
  println("mean_convergence_time")
  println("$mean_convergence_time")

  plot_trajectories(xds, xc, hd, hc, T), plot(1:Hd, dev₁_trajs[1], lab="\$\\Delta x_1\$", lw=2)
end
end

# ╔═╡ afa272df-acd5-4df4-9f96-b000869742e9
md"""
## Combined Edge and Cloud Input
"""

# ╔═╡ a1bba739-9095-4708-92fb-ffef446dd315
let
begin
  z₀ = [fill(x0, size(sys.A, 1) * (n + 1)); u0; u0]
  sys_d, K_d = synthesize(sys, hd, Dcₑ, Dcₔ, α)

  xds, uds = unzip([simulate(sys_d, Hd, z₀, K_uncertain, K_d, σₑ, σₔ, μ) for _ in 1:N])
  xds = [[xd_i[1:2] for xd_i in xd] for xd in xds]
  uds = [[ud_i[1] for ud_i in ud] for ud in uds]
end

begin
  #rel_dev_trajs = [[(1 .- (xds[i][k] ./ xc[1+(k-1)*floor(Int,hd/hc)])) for k in 1:Hd] for i in 1:N]
  dev_trajs = [[(xc[1+(k-1)*floor(Int, hd / hc)] .- xds[i][k]) for k in 1:Hd] for i in 1:N]
  dev₁_trajs = [[dev[1] for dev in dev_traj] for dev_traj in dev_trajs]
  uds₁ = [[udi[1] for udi in ud] for ud in uds]

  max_max_x₁ = maximum(xd -> maximum(xd -> abs(xd[1]), xd), xds)

  mean_mean_dev₁ = mean([mean(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])
  max_max_dev₁ = maximum([maximum(abs.(dev₁_traj)) for dev₁_traj in dev₁_trajs])

  mean_rmse_dev₁ = mean([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  std_rmse_dev₁ = std([std(dev₁_traj) for dev₁_traj in dev₁_trajs])
  max_rmse_dev₁ = maximum([std(dev₁_traj) for dev₁_traj in dev₁_trajs])

  mean_mean_u₁ = mean([mean(abs.(ud₁)) for ud₁ in uds₁])
  std_mean_u₁ = std([mean(abs.(ud₁)) for ud₁ in uds₁])
  max_max_u₁ = maximum([maximum(abs.(ud₁)) for ud₁ in uds₁])

  convergences = [first_convergence([0.0, 0.0], xds[i], threshold=1e-2) for i in 1:N]
  mean_convergence_time = mean(filter(!isnothing, convergences)) * hd

  max_xc₁ = maximum(xc -> abs(xc[1]), xc)
  max_uc₁ = maximum(u -> abs(u[1]), uc)
end

begin
  println("mean_rmse_dev₁	std_rmse_dev₁	mean_mean_u₁	std_mean_u₁")
  println("$mean_rmse_dev₁	$std_rmse_dev₁	$mean_mean_u₁	$std_mean_u₁")
  println("")
  println("max_max_x₁	max_xc₁	max_max_u₁	max_uc₁")
  println("$max_max_x₁	$max_xc₁	$max_max_u₁	$max_uc₁")
  println("")
  println("mean_convergence_time")
  println("$mean_convergence_time")

  plot_trajectories(xds, xc, hd, hc, T), plot(1:Hd, dev₁_trajs[1], lab="\$\\Delta x_1\$", lw=2)
end
end

# ╔═╡ c147b4b5-f8a2-41cf-81bb-36b91d8bff02
md"""
# Miscellaneous
"""

# ╔═╡ ad0f8ae3-ffcd-416f-80a2-7a80483370a0
#println("Trajectories that failed to converge:\n", [i for i in 1:N if isnothing(first_convergence([0.0, 0.0], xds[i], threshold=1e-2))])

# ╔═╡ Cell order:
# ╟─78de028c-1eb1-4120-a3c5-455754ffa9fb
# ╟─6d3ba880-d5f7-4cdf-83db-fd321031f867
# ╟─08bed08d-7ce8-4f07-b055-c1e5b9e33ade
# ╟─cb10d46e-378a-4726-90ea-e1f4b4be796f
# ╟─3cd33b19-87de-4280-8a5f-73f071c08c4b
# ╟─ca5f1534-4201-4dab-a917-9d0d15b0f13d
# ╟─f73455e6-e275-4b02-93fe-cb311a60fdc4
# ╟─565778df-82b2-448d-a6bd-626d8fe7fc31
# ╟─3635e9ac-8429-49bf-93dd-96eec97b807f
# ╟─9e2f99ea-b711-492f-915d-4981325856e7
# ╟─ac159d66-9878-4aca-aa88-0dd576481010
# ╟─24abfc73-8ff4-4a2b-b377-78c134fea835
# ╟─8f06880f-0046-4569-a790-19b03988f36a
# ╟─f421fe0b-b121-4c74-9ea2-b45c328759b6
# ╟─20d653a3-9489-4671-9432-ddad7e83bc8a
# ╟─0d3bfeca-b728-466a-a8c9-4b4187072d55
# ╟─2d2f7f38-fe61-44d4-b3a4-0bd584da7c1a
# ╟─e9b3e724-0703-4aec-a8be-1af5da0a9ebe
# ╟─cb49ecee-ff53-4fa2-98be-e752885de5e2
# ╟─51add47f-1c67-4ac3-a3a2-d5d1f870c1a7
# ╟─672d6769-0077-4c3a-9991-7f70b72edac4
# ╟─c2cca18a-1434-464f-b774-06f217d00271
# ╟─86014118-222d-4cbe-b4ce-b90c6e6d9078
# ╟─68b22a80-d58b-47a5-ad16-83591913839b
# ╟─afa272df-acd5-4df4-9f96-b000869742e9
# ╟─a1bba739-9095-4708-92fb-ffef446dd315
# ╟─c147b4b5-f8a2-41cf-81bb-36b91d8bff02
# ╟─ad0f8ae3-ffcd-416f-80a2-7a80483370a0
