module NoisyReach

export reach, get_error_bound, max_diam
include("reachability.jl")

export benchmarks, c2d
include("models.jl")

export unzip, int_expAs_B, normal_sample, simulate, synthesize, K_uncertain, plot_trajectories
include("helpers.jl")

end
