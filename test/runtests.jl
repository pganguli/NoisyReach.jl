using NoisyReach

using Test
using ControlSystemsBase
using LinearAlgebra
using ReachabilityAnalysis

# Neural network latency
const latency = 0.02
# Error of each individual dimension
const errors = [0.27, 0.27]

@testset "NoisyReach.jl" begin
  # Linearized bicycle model of the F1/10 race car
  sys = let
    v = 6.5
    L = 0.3302
    d = 1.5
    A = [0 v; 0 0]
    B = [0; v / L]
    C = [1 0]
    D = 0

    ss(A, B, C, D)
  end

  A = c2d(sys, latency).A
  B = c2d(sys, latency).B
  K = lqr(ControlSystemsBase.Discrete, A, B, I, I)

  @test size(A) == (2, 2)
  @test size(B) == (2, 1)
  @test size(K) == (1, 2)

  # Perception error bound zonotope
  E = Zonotope(zeros(Float64, 2), Diagonal(errors))

  # Closed-loop dynamics
  Φ = A - B * K

  # Overall error bound zonotope
  W = get_error_bound(B, K, E)
  @test W.center == [0.0, 0.0]
  @test size(W.generators) == (2, 1)
  @test isequivalent(W, -B * K * E)

  x0center = 10.0
  x0size = 1.0
  x0 = Zonotope(x0center * ones(2), x0size * I(2))

  r = reach(Φ, x0, W, 100)
  @time r = reach(Φ, x0, W, 100)
  @info max_diam(r)
  @test max_diam(r) > 0

  # Testing reachability analysis with the error bound zonotope W defined as a function
  W_func = (x) -> let
    # Get the vertices of the zonotope; `stack` combines the list of vectors into a matrix.
    vertices_matrix = vertices_list(x) |> stack
    # Find the maximum value for each state dimension, after taking the absolute values.
    max_states = maximum(abs.(vertices_matrix), dims=2) |> vec
    # Calculate the new error bound zonotope
    get_error_bound(B, K, Zonotope(zeros(Float64, 2), Diagonal(errors .* max_states)))
  end
  # For some reason this is really slow; I'm setting max_order to 10 for now.
  @time r = reach(Φ, x0, W_func, 100, max_order=10)
  @info max_diam(r)
  @test max_diam(r) > 0
end
