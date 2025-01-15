using QuadGK: quadgk
using Distributions: Normal
using LinearAlgebra: I, ℯ
using ControlSystemsBase: ss, lqr, Continuous, Discrete, StateSpace
using Plots

unzip(a) = (getfield.(a, x) for x in fieldnames(eltype(a)))

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
Also note that, ``λ \\in [-1, 1]``.
"""
function normal_sample(σ::Real, μ::Real)
  𝒩 = Normal(μ, σ)
  λ = clamp(rand(𝒩), -1, 1)
  return λ
end

"""
	simulate(sys, K, H, x₀)

Return `x, u`, where ``x[k+1] = sys.A \\cdot x[k] + sys.B \\cdot u[k]`` and ``u[k+1] = K \\cdot x[k+1]`` for `H` steps.
Also note that, ``x[0] = x₀``.
"""
function simulate(sys::StateSpace, K::AbstractMatrix, H::Integer, x₀::AbstractVector)
  x = Vector{Vector{Float64}}(undef, H + 1)
  u = Vector{Vector{Float64}}(undef, H)

  x[1] = x₀
  for k in 1:H
    uₖ = K * x[k]; u[k] = isa(uₖ, AbstractVector) ? uₖ : [uₖ]
    x[k+1] = sys.A * x[k] + sys.B * u[k]
  end

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
	synthesize(sys, h)

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
	synthesize(sys, h)

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
	K_uncertain(K, σ, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated by `σ` centered around `μ`.
"""
function K_uncertain(K::AbstractMatrix, σ::Real, μ::Real)
  λ = [normal_sample(σ, μ) for _ in 1:size(K, 2)-1]
  K_ = [(K[1, 1:end-1] .* (1 .+ λ)); K[1, end]]'
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
    [(K[1, 1:end-1] .* (1 .+ λ₁)); K[1, end]]';
    [(K[2, 1:end-1] .* (1 .+ λ₂)); K[2, end]]'
  ]
  return K_
end

function plot_trajectories(traj::AbstractVector, traj_ideal::AbstractVector, xlim::Real, ylim::Real; fname::String="", xlabel::String="\$x_1\$", ylabel="\$x_2\$", title::String="", fontsize::Real=1.8)
  Plots.scalefontsizes(fontsize)

  traj_plot = plot(xlabel=xlabel, ylabel=ylabel, title=title)

  for tr in traj
      plot!([pt[1] for pt in tr], [pt[2] for pt in tr],
            xlim=(0, xlim), ylim=(0, ylim), label="", linecolor=:lightgray, linewidth=1)
  end
  plot!([pt[1] for pt in traj_ideal], [pt[2] for pt in traj_ideal],
        xlim=(0, xlim), ylim=(0, ylim), label="", linecolor=:black, linewidth=2, marker=:circle, markercolor=:red, markersize=3)

  fname != "" && savefig(traj_plot, fname)

  return traj_plot
end
