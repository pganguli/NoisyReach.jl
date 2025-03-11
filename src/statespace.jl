using QuadGK: quadgk
using LinearAlgebra: I, ℯ, kron
using ControlSystemsBase: ss, lqr, Continuous, Discrete, StateSpace

"""
	int_expAs_B(A, B, lo, hi)

Compute ``\\int_{\\text{lo}}^{\\text{hi}} e^{As} \\cdot B \\, ds``.
"""
function int_expAs_B(A::AbstractMatrix, B::AbstractMatrix, lo::Real, hi::Real)
  return first(quadgk(s -> exp(A * s) * B, lo, hi))
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
    uₖ = let
      tmp = -K(args...) * x[k]
      isa(tmp, AbstractVector) ? tmp : [tmp]
    end
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

Synthesize discrete-time state space model `_sys` and corresponding controller `K` for continuous-time state space model `sys`,
for sampling period `h`.
NOTE: Usual discretization.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real)
  _sys = let
    ϕ = ℯ^(h * sys.A)
    Γ = int_expAs_B(sys.A, sys.B, 0, h)
    ss(ϕ, Γ, sys.C, sys.D, h)
  end
  K = lqr(Discrete, _sys.A, _sys.B, I, I)
  return _sys, K
end

"""
	synthesize(sys, h, Dc)

Synthesize discrete-time state space model `_sys` and corresponding controller `K` for continuous-time state space model `sys`,
for sampling period `h` and input-delay `Dc`.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc::Real)
  n = floor(Int, Dc / h)
  _Dc = Dc - n * h
  _sys = let
    ϕ = ℯ^(h * sys.A)
    Γ₁ = int_expAs_B(sys.A, sys.B, h - _Dc, h)
    Γ₂ = int_expAs_B(sys.A, sys.B, 0, h - _Dc)
    ϕ_aug = let
      I_block = let
        tmp = kron(I(n), I(size(ϕ, 2)))
        size(tmp, 2) != 0 ? tmp : zeros(0, size(ϕ, 2))
      end
      Z_row = zeros(1, size(ϕ, 2) * (n + 1) + size(Γ₁, 2))
      Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
      ϕ_row = [ϕ zeros(size(ϕ, 1), size(Z_row, 2) - size(Gamma_col₁, 2) - size(ϕ, 2))]
      Z_col = zeros(size(ϕ, 1) * (n), size(Z_row, 2) - size(Gamma_col₁, 2) - size(I_block, 2))
      vcat(hcat(vcat(ϕ_row, hcat(I_block, Z_col)), Gamma_col₁), Z_row)
    end
    Γ_aug = [Γ₂; zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I]
    C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
    D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
    ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
  end
  K = lqr(Discrete, _sys.A, _sys.B, I, I)
  return _sys, K
end

"""
	synthesize(sys, h, Dc₁, Dc₂, α)

Synthesize discrete-time state space model `_sys` and corresponding controller `K` for continuous-time state space model `sys`,
for sampling period `h` and input-delays `Dc₁` and `Dc₂`, merged in `α` ratio.
"""
function synthesize(sys::StateSpace{Continuous}, h::Real, Dc₁::Real, Dc₂::Real, α::Real)
  n = floor(Int, Dc₂ / h)
  _Dc₂ = Dc₂ - n * h
  ue_before_uc = (Dc₂ - n * h) > Dc₁
  if ue_before_uc
    # NOTE: Split computing delays assuming `u_c[k-n]` arrives after `u_e[k]`.
    _sys = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, h - _Dc₂, h - _Dc₂ + Dc₁)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0, _Dc₂ - Dc₁)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0, h - _Dc₂)
      ϕ_aug = let
        I_block = let
          tmp = kron(I(n), I(size(ϕ, 2)))
          size(tmp, 2) != 0 ? tmp : zeros(0, size(ϕ, 2))
        end
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + 2 * size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        Gamma_col₂ = [Γ₁ * (1 - α); zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(Z_row, 2) - size(Gamma_col₁, 2) - size(Gamma_col₂, 2) - size(ϕ, 2))]
        Z_col = zeros(size(ϕ, 1) * (n), size(Z_row, 2) - size(Gamma_col₁, 2) - size(Gamma_col₂, 2) - size(I_block, 2))
        vcat(hcat(vcat(ϕ_row, hcat(I_block, Z_col)), Gamma_col₁, Gamma_col₂), Z_row, Z_row)
      end
      Γ_aug = [(Γ₃*α+Γ₂) (Γ₃*(1-α)); zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)) zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)); I zeros(1, size(Γ₃, 2)); zeros(1, size(Γ₂, 2)) I]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, _sys.A, _sys.B, I, I)
    return _sys, K
  else
    # NOTE: Split computing delays assuming `u_c[k-n]` arrives before `u_e[k]`.
    _sys = let
      ϕ = ℯ^(h * sys.A)
      Γ₁ = int_expAs_B(sys.A, sys.B, h - Dc₁, h - Dc₁ + _Dc₂)
      Γ₂ = int_expAs_B(sys.A, sys.B, 0, Dc₁ - _Dc₂)
      Γ₃ = int_expAs_B(sys.A, sys.B, 0, h - Dc₁)
      ϕ_aug = let
        I_block = let
          tmp = kron(I(n), I(size(ϕ, 2)))
          size(tmp, 2) != 0 ? tmp : zeros(0, size(ϕ, 2))
        end
        Z_row = zeros(1, size(ϕ, 2) * (n + 1) + 2 * size(Γ₁, 2))
        Gamma_col₁ = [Γ₁; zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        Gamma_col₂ = [Γ₁ * (1 - α); zeros(size(Γ₁, 1) * (n), size(Γ₁, 2))]
        ϕ_row = [ϕ zeros(size(ϕ, 1), size(Z_row, 2) - size(Gamma_col₁, 2) - size(Gamma_col₂, 2) - size(ϕ, 2))]
        Z_col = zeros(size(ϕ, 1) * (n), size(Z_row, 2) - size(Gamma_col₁, 2) - size(Gamma_col₂, 2) - size(I_block, 2))
        vcat(hcat(vcat(ϕ_row, hcat(I_block, Z_col)), Gamma_col₁, Gamma_col₂), Z_row, Z_row)
      end
      Γ_aug = [(Γ₃*α) (Γ₃*(1-α)+Γ₂); zeros(size(Γ₃, 1) * (n), size(Γ₃, 2)) zeros(size(Γ₂, 1) * (n), size(Γ₂, 2)); I zeros(1, size(Γ₂, 2)); zeros(1, size(Γ₃, 2)) I]
      C_aug = [sys.C zeros(size(sys.C, 1), size(ϕ_aug, 2) - size(sys.C, 2))]
      D_aug = [sys.D zeros(size(sys.D, 1), size(Γ_aug, 2) - size(sys.D, 2))]
      ss(ϕ_aug, Γ_aug, C_aug, D_aug, h)
    end
    K = lqr(Discrete, _sys.A, _sys.B, I, I)
    return _sys, K
  end
end
