using QuadGK: quadgk
using LinearAlgebra: I, ℯ, kron
using ControlSystemsBase: ss, lqr, Continuous, Discrete, StateSpace

"""
	int_expAs_B(A, B, lo, hi)

Compute ``\\int_{\\text{lo}}^{\\text{hi}} e^{As} \\cdot B \\, ds``.
"""
function int_expAs_B(A::AbstractMatrix, B::AbstractMatrix, lo::Real, hi::Real)
  integrand = s -> exp(A * s) * B
  integral, _ = quadgk(integrand, lo, hi)
  return integral
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
  for k in 1:H
    uₖ = -K(args...) * x[k]
    uₖ = isa(uₖ, AbstractVector) ? uₖ : [uₖ]
    if any(ismissing.(uₖ))
      u[k] = k == 1 ? zeros(size(uₖ)) : u[k-1]
    else
      u[k] = uₖ
    end
    x[k+1] = sys.A * x[k] + sys.B * u[k]
  end

  return x, u
end

"""
	synthesize(sys, h)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h`.
NOTE: Usual discretization.
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
NOTE: [Case-3] Usual discretization with computation delay.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₀ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₁)
    Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc₁, h)
    ϕ_aug = [ϕ Γ₁; zeros(1, size(ϕ, 2) + size(Γ₁, 2))]
    Γ_aug = [Γ₀; I]
    C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
    D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc₁, Dc₂)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
NOTE: [Case-0] Split computing delays assuming `Dc₁`, `Dc₂` < `h`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₃ = int_expAs_B(sys.A, sys.B, h - Dc₂, h - Dc₂ + Dc₁)
    Γ₁ = int_expAs_B(sys.A, sys.B, 0.0, Dc₂ - Dc₁)
    Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₂)
    ϕ_aug = [ϕ Γ₃; zeros(1, size(ϕ, 2) + size(Γ₃, 2))]
    Γ_aug = [Γ₁ Γ₂; zeros(1, size(Γ₁, 2)) I]
    C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
    D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delay `Dc`.
NOTE: [Case-4] Usual discretization with computation delay > `n*h`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₂::Real, n::Integer)
  sys_ = let
    ϕ = ℯ^(h * sys.A)
    Γ₁ = int_expAs_B(sys.A, sys.B, (n + 1) * h - Dc₂, Dc₂ - n * h)
    Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, (n + 1) * h - Dc₂)
    ϕ_aug = let
      I_block = kron(I(n), I(size(ϕ, 2)))
      Z_col = zeros(size(ϕ, 1) * (n + 1), size(ϕ, 2))
      Z_row = zeros(1, size(ϕ, 2) * (n + 1) + size(Γ₁, 2))
      Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
      ϕ_row = [ϕ zeros(size(ϕ, 1), size(ϕ, 2) * (n - 1))]
      vcat(hcat(vcat(ϕ_row, I_block), Z_col, Gamma_col₁), Z_row)
    end
    Γ_aug = [Γ₂; zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I]
    C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
    D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, sys_.A, sys_.B, I, I)
  return sys_, K
end

"""
	synthesize(sys, h, Dc₁, Dc₂, n)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
NOTE: [Case-1] Split computing delays assuming `u_c[k-n]` arrives before `u_e[k]`.
NOTE: [Case-2] Split computing delays assuming `u_e[k]` arrives before `u_c[k-n]`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real, n::Integer; ue_before_uc::Bool)
  # NOTE: [Case-2] Split computing delays assuming `u_e[k]` arrives before `u_c[k-n]`.
  if ue_before_uc
    sys_ = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, (n + 1) * h - Dc₂, (n + 1) * h - Dc₂ + Dc₁)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, Dc₂ - Dc₁ - n * h)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0.0, (n + 1) * h - Dc₂)
      ϕ_aug = let
        I_block = kron(I(n), I(size(ϕ, 2)))
        Z_col = zeros(size(ϕ, 1) * (n + 1), size(ϕ, 2))
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(ϕ, 2) * (n - 1))]
        vcat(hcat(vcat(ϕ_row, I_block), Z_col, Gamma_col₁), Z_row)
      end
      Γ_aug = [Γ₃ Γ₂; zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)) zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I zeros(1, size(Γ₂, 2))]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, sys_.A, sys_.B, I, I)
    return sys_, K
    # NOTE: [Case-1] Split computing delays assuming `u_c[k-n]` arrives before `u_e[k]`.
  else
    sys_ = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc₁, Dc₂ - Dc₁ - (n - 1) * h)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, n * h - Dc₂ + Dc₁)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₁)
      ϕ_aug = let
        I_block = kron(I(n), I(size(ϕ, 2)))
        Z_col = zeros(size(ϕ, 1) * (n + 1), size(ϕ, 2))
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(ϕ, 2) * (n - 1))]
        vcat(hcat(vcat(ϕ_row, I_block), Z_col, Gamma_col₁), Z_row)
      end
      Γ_aug = [Γ₃ Γ₂; zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)) zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I zeros(1, size(Γ₂, 2))]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, sys_.A, sys_.B, I, I)
    return sys_, K
  end
end

"""
	synthesize(sys, h, Dc₁, Dc₂, n)

Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
NOTE: [Case-5] Split computing delays assuming `u_e[k]` arrives before `u_c[k-n]`, merge with α ratio.
NOTE: [Case-6] Split computing delays assuming `u_c[k-n]` arrives before `u_e[k]`, merge with α ratio.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real, n::Integer, α::Real; ue_before_uc::Bool)
  # NOTE: [Case-5] Split computing delays assuming `u_c[k-n]` arrives after `u_e[k]`.
  if ue_before_uc
    sys_ = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, (n + 1) * h - Dc₂, (n + 1) * h - Dc₂ + Dc₁)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, Dc₂ - Dc₁ - n * h)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0.0, (n + 1) * h - Dc₂)
      ϕ_aug = [ϕ Γ₁ Γ₁*(1-α); 0 0 0 0; 0 0 0 0]
      ϕ_aug = let
        I_block = kron(I(n), I(size(ϕ, 2)))
        Z_col = zeros(size(ϕ, 1) * (n + 1), size(ϕ, 2))
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + 2 * size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        Gamma_col₂ = [Γ₁ * (1 - α); zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(ϕ, 2) * (n - 1))]
        vcat(hcat(vcat(ϕ_row, I_block), Z_col, Gamma_col₁, Gamma_col₂), Z_row, Z_row)
      end
      Γ_aug = [(Γ₃*α+Γ₂) (Γ₃*(1-α)); zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)) zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)); I zeros(1, size(Γ₃, 2)); zeros(1, size(Γ₂, 2)) I]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, sys_.A, sys_.B, I, I)
    return sys_, K
    # NOTE: [Case-6] Split computing delays assuming `u_c[k-n]` arrives before `u_e[k]`.
  else
    sys_ = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc₁, Dc₂ - Dc₁ - (n - 1) * h)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0.0, n * h - Dc₂ + Dc₁)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc₁)
      ϕ_aug = let
        I_block = kron(I(n), I(size(ϕ, 2)))
        Z_col = zeros(size(ϕ, 1) * (n + 1), size(ϕ, 2))
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + 2 * size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        Gamma_col₂ = [Γ₁ * (1 - α); zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(ϕ, 2) * (n - 1))]
        vcat(hcat(vcat(ϕ_row, I_block), Z_col, Gamma_col₁, Gamma_col₂), Z_row, Z_row)
      end
      Γ_aug = [(Γ₃*α) (Γ₃*(1-α)+Γ₂); zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)) zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I zeros(1, size(Γ₂, 2)); zeros(1, size(Γ₃, 2)) I]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, sys_.A, sys_.B, I, I)
    return sys_, K
  end
end

#"""
#	synthesize(sys, h, Dc₁, Dc₂, n)
#
#Synthesize discrete-time state space model `sys_` and corresponding controller `K` for conntinuous-time state space model `sys`,
#for sampling period `h` and input-delays `Dc₁` and `Dc₂`.
#NOTE: [Case-7] Split computing delays assuming `u_e[k]` merge with `u_c[k-n]` at α ratio.
#"""
#function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real, n::Integer, α::Real)
#  sys_ = let
#    ϕ = ℯ^(h * sys.A)
#    Γ₀ = int_expAs_B(sys.A, sys.B, 0.0, h - Dc)
#    Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc, h)
#    ϕ_aug = [ϕ Γ₁; 0 0 0]
#    Γ_aug = [Γ₀; I]
#    C_aug = [sys.C 0]
#    ss(ϕ_aug, Γ_aug, C_aug, sys.D, h)
#  end
#  K = lqr(Discrete, sys_.A, sys_.B, I, I)
#  return sys_, K
#end
