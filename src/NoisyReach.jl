module NoisyReach

export reach, get_error_bound, max_diam
include("reachability.jl")

export benchmarks, c2d
include("models.jl")

export unzip, int_expAs_B, normal_sample, iid_sample, seq_sample, simulate, synthesize, K_certain, K_uncertain, plot_trajectories, plot_trajectory, has_converged, first_convergence
include("helpers.jl")

end
