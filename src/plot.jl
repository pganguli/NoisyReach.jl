using Statistics: mean
using Plots: plot, plot!, xlims!, ylims!, savefig

"""
    plot_trajectories(trajs, traj_ideal; xlim=nothing, ylim=nothing, fname="", xlabel="\$x_1\$", ylabel="\$x_2\$", title="")

Plot the trajectories with optional axis limits and save to file if `fname` is provided.
"""
function plot_trajectories(
  trajs::Vector{Vector{Vector{Float64}}}, traj_ideal::Vector{Vector{Float64}}, hd::Real, hc::Real, T::Real;
  xlim::Union{Real,Nothing}=nothing, ylim::Union{Real,Nothing}=nothing,
  fname::String="", title::String="", plot_y=false)
  traj_plot = plot(xlabel="\$Time [s]\$", title=title, fmt=:svg)

  data_to_plot = hcat(traj_ideal...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hc:T, data_to_plot, lab=["x₁(t)" "x₂(t)"], linecolor=[:red :blue], alpha=0.6, fmt=:svg)

  max_traj = [maximum(traj[i] for traj in trajs) for i in 1:length(trajs[1])]
  min_traj = [minimum(traj[i] for traj in trajs) for i in 1:length(trajs[1])]
  avg_traj = [mean(traj[i] for traj in trajs) for i in 1:length(trajs[1])]
  data_to_plot = hcat(max_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab=["max. x₁[k]" "max. x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, linestyle=:dash, fmt=:svg)
  data_to_plot = hcat(min_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab=["min. x₁[k]" "min. x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, linestyle=:dash, fmt=:svg)
  data_to_plot = hcat(avg_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab=["avg. x₁[k]" "avg. x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, linestyle=:dot, fmt=:svg)

  if !isnothing(xlim) && !isnothing(ylim)
    xlims!(0, xlim)
    ylims!(0, ylim)
  end

  if !isempty(fname)
    savefig(traj_plot, fname)
  end

  return traj_plot
end

#"""
#    plot_trajectories(traj, traj_ideal, hc, hd, T; xlim=nothing, ylim=nothing, fname="", xlabel="\$x_1\$", ylabel="\$x_2\$", title="")
#
#Plot the trajectories with optional axis limits and save to file if `fname` is provided.
#"""
#function plot_trajectories(
#  trajs::Vector{Vector{Float64}}, traj_ideal::Vector{Vector{Float64}}, hd::Real, hc::Real, T::Real;
#  xlim::Union{Real,Nothing}=nothing, ylim::Union{Real,Nothing}=nothing,
#  fname::String="", title::String="")
#  traj_plot = plot(xlabel="\$Time [s]\$", title=title, fmt=:svg)
#  traj_plot = plot!(0:hc:T, hcat(traj_ideal...)', lab=["x₁(t)" "x₂(t)"], linecolor=[:red :blue], alpha=0.6, fmt=:svg)
#  traj_plot = plot!(0:hd:T, hcat(traj...)', lab=["x₁[k]" "x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, marker=:circle, markersize=2, fmt=:svg)
#
#  if !isnothing(xlim) && !isnothing(ylim)
#    xlims!(0, xlim)
#    ylims!(0, ylim)
#  end
#
#  if !isempty(fname)
#    savefig(traj_plot, fname)
#  end
#
#  return traj_plot
#end
