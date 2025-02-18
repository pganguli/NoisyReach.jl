using QuadGK: quadgk
using Distributions: Normal
using LinearAlgebra: I, ℯ
using ControlSystemsBase: ss, lqr, Continuous, Discrete, StateSpace
using Plots

unzip(a) = (getfield.(a, x) for x in fieldnames(eltype(a)))
global COUNTER = 1

"""
	int_expAs_B(A, B, lo, hi)

Compute ``\\int_{\\text{lo}}^{\\text{hi}} e^{As} \\cdot B \\, ds``.
"""
function int_expAs_B(A::AbstractMatrix, B::AbstractMatrix, lo::Real, hi::Real)
  integrand = s -> exp(A * s) * B
  result, _ = quadgk(integrand, lo, hi)
  return result
end

"""
	normal_sample(σ, μ)

Returns `λ`, where ``λ \\sim 𝒩(μ, σ^2)``.
"""
function normal_sample(σ::Real, μ::Real)
  𝒩 = Normal(μ, σ)
  λ = rand(𝒩)
  @debug "λ=$λ"
  return λ
end

"""
	iid_sample(experiment)

Returns `λ`, where ``λ \\sim i.i.d(experimental data)``
"""
function iid_sample(list::AbstractVector)
  λ = rand(list)
  @debug "λ=$λ"
  return λ
end

"""
	seq_sample(experiment)

Returns `λ`, where ``λ \\sim seq(experimental data)``
"""
function seq_sample(list::AbstractVector)
  global COUNTER
  λ = list[COUNTER % size(list,1)]
  COUNTER += 1
  @debug "λ=$λ\tCOUNTER=$COUNTER"
  return λ
end

"""
	simulate(sys, K, H, x₀)

Return `x, u`, where ``x[k+1] = sys.A \\cdot x[k] + sys.B \\cdot u[k]`` and ``u[k+1] = K \\cdot x[k+1]`` for `H` steps.
Also note that, ``x[0] = x₀``.
"""
function simulate(sys::StateSpace, H::Integer, x₀::AbstractVector, K, args...)
  x = Vector{Vector{Float64}}(undef, H + 1)
  u = Vector{Vector{Float64}}(undef, H)

  x[1] = x₀
  @debug "---SIMULATION BEGIN---"
  for k in 1:H
    uₖ = -K(args...) * x[k]
    uₖ = isa(uₖ, AbstractVector) ? uₖ : [uₖ]
    if any(ismissing.(uₖ))
      u[k] = k == 1 ? zeros(size(uₖ)) : u[k-1]
    else
      u[k] = uₖ
    end
    @debug "Using u[$k]=$(u[k])"
    x[k+1] = sys.A * x[k] + sys.B * u[k]
  end
  @debug "---SIMULATION END---"

  return x, u
end

"""
	synthesize(sys, h)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ = int_expAs_B(sys.A, sys.B, 0.0, h)
    ss(ϕ, Γ, sys.C, sys.D, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delay `Dc`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc::Real)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₀ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc)
    Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc, h)
    ϕ_aug = [ϕ Γ₁; 0 0 0]
    Γ_aug = [Γ₀; I]
    C_aug = [sys.C 0]
    D_aug = sys.D
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc₁, Dc₂)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₃ = int_expAs_B(sys.A, sys.B, h - Dc₂, h - Dc₂ + Dc₁)
    Γ₁ = int_expAs_B(sys.A, sys.B, 0.0, Dc₂ - Dc₁)
    Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₂)
    ϕ_aug = [ϕ Γ₃; 0 0 0]
    Γ_aug = [Γ₁ Γ₂; 0 I]
    C_aug = [sys.C 0]
    D_aug = [sys.D 0]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc₁, Dc₂, n)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real, n::Integer)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc₁, Dc₂ - Dc₁ - (n - 1) * h)
    Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, n * h - Dc₂ + Dc₁)
    Γ₃ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₁)
    ϕ_aug = [ϕ Γ₁; 0 0 0]
    Γ_aug = [Γ₃ Γ₂; 0 I]
    C_aug = [sys.C 0]
    D_aug = [sys.D 0]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	K_certain(K)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is not deviated.
"""
function K_certain(K::AbstractMatrix)
  return K
end

"""
	K_uncertain(K, σ, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated as per i.i.d assumption on experimental data.
"""
function K_uncertain(K::AbstractMatrix, list::AbstractVector)
  λ = [seq_sample(list) for _ in 1:size(K, 2)-1]
  K_ = [(K[1, 1:end-1] .* (1 .- λ)); K[1, end]]'
  @debug "λ=$λ"
  return K_
end

"""
	K_uncertain(K, σ, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated by `σ` centered around `μ`.
"""
function K_uncertain(K::AbstractMatrix, σ::Real, μ::Real)
  λ = [normal_sample(σ, μ) for _ in 1:size(K, 2)-1]
  K_ = [(K[1, 1:end-1] .* (1 .- λ)); K[1, end]]'
  @debug "λ=$λ"
  return K_
end

"""
	K_uncertain(K, σ₁, σ₂, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated by `σ₁, σ₂` centered around `μ`.
"""
function K_uncertain(K::AbstractMatrix, σ₁::Real, σ₂::Real, μ::Real)
  λ₁ = [normal_sample(σ₁, μ) for _ in 1:size(K, 2)-1]
  λ₂ = [normal_sample(σ₂, μ) for _ in 1:size(K, 2)-1]
  K_ = [
    [(K[1, 1:end-1] .* (1 .- λ₁)); K[1, end]]';
    [(K[2, 1:end-1] .* (1 .- λ₂)); K[2, end]]'
  ]
  @debug "λ₁=$λ₁"
  @debug "λ₂=$λ₂"
  return K_
end

"""
	plot_trajectories(traj, traj_ideal, xlim, ylim; fname, xlabel, ylabel, title)

Plot the trajectories.
"""
function plot_trajectories(traj::AbstractVector, traj_ideal::AbstractVector; xlim::Real=0, ylim::Real=0, fname::String="", xlabel::String="\$x_1\$", ylabel="\$x_2\$", title::String="")
  traj_plot = plot(xlabel=xlabel, ylabel=ylabel, title=title)

  for tr in traj
    if (xlim != 0 || ylim != 0)
      plot!([pt[1] for pt in tr], [pt[2] for pt in tr],
        xlim=(0, xlim), ylim=(0, ylim), label="", linecolor=:lightgray, linewidth=1)
    else
      plot!([pt[1] for pt in tr], [pt[2] for pt in tr],
        label="", linecolor=:lightgray, linewidth=1)
    end
  end
  if (xlim != 0 || ylim != 0)
    plot!([pt[1] for pt in traj_ideal], [pt[2] for pt in traj_ideal],
      xlim=(0, xlim), ylim=(0, ylim), label="", linecolor=:black, linewidth=2, marker=:circle, markercolor=:red, markersize=3)
  else
    plot!([pt[1] for pt in traj_ideal], [pt[2] for pt in traj_ideal],
      label="", linecolor=:black, linewidth=2, marker=:circle, markercolor=:red, markersize=3)
  end

  fname != "" && savefig(traj_plot, fname)

  return traj_plot
end

"""
	plot_trajectory(xds, hd, T, i)

Plot the i-th trajectory.
"""
function plot_trajectory(xds::AbstractVector, hd::Real, T::Real, i::Integer)
  plot!(0:hd:T, hcat(xds[i]...)', lab=["x₁[k]" "x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, marker=:circle, markersize=2)
end

function has_converged(x_ref::Vector{<:Real}, x::Vector{<:Real}; threshold::Real=0.01)
  return all(-(threshold .- x_ref) .<= x .<= (threshold .+ x_ref))
end

function first_convergence(x_ref::AbstractVector, xs::AbstractVector; threshold::Real=0.01)
  return findfirst(x -> has_converged(x_ref, x, threshold=threshold), xs)
end
