using Distributions: Normal

unzip(a) = (getfield.(a, x) for x in fieldnames(eltype(a)))
global COUNTER = 0

"""
	normal_sample(σ, μ)

Returns `λ`, where ``λ \\sim 𝒩(μ, σ^2)``.
"""
function normal_sample(σ::Real, μ::Real)
  𝒩 = Normal(μ, σ)
  λ = rand(𝒩)
  return λ
end

"""
	iid_sample(experiment)

Returns `λ`, where ``λ \\sim i.i.d(experimental data)``
"""
function iid_sample(list::AbstractVector)
  λ = rand(list)
  return λ
end

"""
	seq_sample(experiment)

Returns `λ`, where ``λ \\sim seq(experimental data)``
"""
function seq_sample(list::AbstractVector)
  global COUNTER
  λ = list[COUNTER % size(list,1) + 1]
  COUNTER += 1
  return λ
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
NOTE: Assume ``x' = (1-λ) \\cdot x``.
"""
function K_uncertain(K::AbstractMatrix, list::AbstractVector)
  λ = [iid_sample(list) for _ in 1:size(K, 2)-1]
  K_ = [(K[1, 1:end-1] .* (1 .- λ)); K[1, end]]'
  return K_
end

"""
	K_uncertain(K, σ, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated by `σ` centered around `μ`.
NOTE: Assume ``x' = (1-λ) \\cdot x``.
"""
function K_uncertain(K::AbstractMatrix, σ::Real, μ::Real)
  λ = [normal_sample(σ, μ) for _ in 1:size(K, 2)-1]
  K_ = [(K[1, 1:end-1] .* (1 .- λ)); K[1, end]]'
  return K_
end

"""
	K_uncertain(K, σ₁, σ₂, μ)

Compute `K'` such that ``u[k+1] = K' \\cdot x[k+1]`` is deviated by `σ₁, σ₂` centered around `μ`.
NOTE: Assume ``x' = (1-λ) \\cdot x``.
"""
function K_uncertain(K::AbstractMatrix, σ₁::Real, σ₂::Real, μ::Real)
  λ₁ = [normal_sample(σ₁, μ) for _ in 1:size(K, 2)-1]
  λ₂ = [normal_sample(σ₂, μ) for _ in 1:size(K, 2)-1]
  K_ = [
    [(K[1, 1:end-1] .* (1 .- λ₁)); K[1, end]]';
    [(K[2, 1:end-1] .* (1 .- λ₂)); K[2, end]]'
  ]
  return K_
end


function has_converged(x_ref::Vector{<:Real}, x::Vector{<:Real}; threshold::Real=0.01)
  return all(-(threshold .- x_ref) .<= x .<= (threshold .+ x_ref))
end

function first_convergence(x_ref::AbstractVector, xs::AbstractVector; threshold::Real=0.01)
  return findfirst(x -> has_converged(x_ref, x, threshold=threshold), xs)
end
