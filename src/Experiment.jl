module Experiment
export matrix_integral, evolve, generate_uncertainty, ideal_evolve

using QuadGK
using Distributions
using ControlSystemsBase
using LinearAlgebra
using ReachabilityAnalysis

function matrix_integral(A::AbstractMatrix, B::AbstractMatrix, lower_limit::Float64, upper_limit::Float64)
  integrand = s -> exp(A * s) * B
  result, _ = quadgk(integrand, lower_limit, upper_limit)
  return result
end

function generate_uncertainty(σ1::Float64, σ2::Float64, μ::Float64)
  dist_λ1 = Normal(μ, σ1)
  dist_λ2 = Normal(μ, σ2)
  λ11 = clamp(rand(dist_λ1), -1, 1)
  λ12 = clamp(rand(dist_λ1), -1, 1)
  λ21 = clamp(rand(dist_λ2), -1, 1)
  λ22 = clamp(rand(dist_λ2), -1, 1)
  return [λ11 λ12 λ21 λ22]
end

function evolve(A::AbstractMatrix, B::AbstractMatrix, K::AbstractMatrix, H::Integer, z0::Vector{Float64}, u1_0::Float64, u2_0::Float64, σ1::Float64, σ2::Float64, μ::Float64)
  u0 = [u1_0; u2_0]
  z = Vector{typeof(z0)}(undef, H + 1)
  u = Vector{typeof(u0)}(undef, H)
  z[1] = z0
  u[1] = u0
  z[2] = A * z[1] + B * u[1]
  z[2][end-length(u2_0)+1:end, :] .= 0
  λ = generate_uncertainty(σ1, σ2, μ)
  λ11 = λ[1]
  λ12 = λ[2]
  λ21 = λ[3]
  λ22 = λ[4]
  K_error = [K[1, :][1]*(1+λ11) K[1, :][2]*(1+λ12) K[1, :][3]; K[2, :][1]*(1+λ21) K[2, :][2]*(1+λ22) K[2, :][3]]
  for k in 2:H
    u[k] = K_error * z[k]
    z[k+1] = A * z[k] + B * u[k]
  end
  return z
end

function ideal_evolve(A::AbstractMatrix, B::AbstractMatrix, K::AbstractMatrix, H::Integer, z0::Vector{Float64}, u1_0::Float64, u2_0::Float64)
  u0 = [u1_0; u2_0]
  z = Vector{typeof(z0)}(undef, H + 1)
  u = Vector{typeof(u0)}(undef, H)
  z[1] = z0
  u[1] = u0
  z[2] = A * z[1] + B * u[1]
  z[2][end-length(u2_0)+1:end, :] .= 0
  for k in 2:H
    u[k] = K * z[k]
    z[k+1] = A * z[k] + B * u[k]
  end
  return z
end
end
