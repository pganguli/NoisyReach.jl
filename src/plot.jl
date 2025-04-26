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
  traj_plot = plot(xlabel="\$Time\\ (s)\$", ylabel="\$x_1\$", title=title, fmt=:svg, guidefont=20, tickfont=16)

  #
  # IDEAL
  #
  data_to_plot = hcat(traj_ideal...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hc:T, data_to_plot, lab="Ideal \$x_1\$", linecolor=:green,   alpha=0.6, linestyle=:solid, marksershape=:none, lw=4, ms=5, mcolor=:green, legendfontsize=16, fmt=:svg)

  max_traj = [maximum(traj[i] for traj in trajs) for i in 1:length(trajs[1])]
  min_traj = [minimum(traj[i] for traj in trajs) for i in 1:length(trajs[1])]
  avg_traj = [mean(traj[i] for traj in trajs) for i in 1:length(trajs[1])]

  #
  # MAX
  #
  data_to_plot = hcat(max_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab="Max. \$x_1\$", linecolor=:darkred,  alpha=0.6, linestyle=:dash,    markershape=:xcross, lw=4, ms=5, mcolor=:darkred, legendfontsize=16, fmt=:svg)

  #
  # MIN
  #
  data_to_plot = hcat(min_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab="Min. \$x_1\$", linecolor=:purple3,     alpha=0.6, linestyle=:dash,    markershape=:cross, lw=4, ms=5, mcolor=:purple3, legendfontsize=16, fmt=:svg)

  #
  # AVG
  #
  data_to_plot = hcat(avg_traj...)'[1:end, plot_y ? (1:2) : 1]
  traj_plot = plot!(0:hd:T, data_to_plot, lab="Avg. \$x_1\$", linecolor=:blue,    alpha=0.6, linestyle=:dashdot, markershape=:circle, lw=4, ms=5, mcolor=:blue, legendfontsize=16, fmt=:svg)

  if !isnothing(xlim)
    xlims!(0, xlim)
  end

  if !isnothing(ylim)
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
