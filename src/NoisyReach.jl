module NoisyReach

export reach, get_error_bound, max_diam
include("reachability.jl")

export benchmarks
include("models.jl")

export int_expAs_B, simulate, synthesize
include("statespace.jl")

export plot_trajectories
include("plot.jl")

export unzip, normal_sample, iid_sample, seq_sample, K_certain, K_uncertain, has_converged, first_convergence
include("helpers.jl")

end
