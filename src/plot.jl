using Plots: plot, plot!, xlims!, ylims!, savefig

"""
    plot_trajectories(trajs, traj_ideal; xlim=nothing, ylim=nothing, fname="", xlabel="\$x_1\$", ylabel="\$x_2\$", title="")

Plot the trajectories with optional axis limits and save to file if `fname` is provided.
"""
function plot_trajectories(trajs::AbstractVector, traj_ideal::AbstractVector;
  xlim::Union{Real,Nothing}=nothing, ylim::Union{Real,Nothing}=nothing,
  fname::String="", title::String="")
  traj_plot = plot(xlabel="\$x_1\$", ylabel="\$x_2\$", title=title)

  for traj in trajs
    plot!([pt[1] for pt in traj], [pt[2] for pt in traj],
      label="", linecolor=:lightgray, linewidth=1)
  end

  plot!([pt[1] for pt in traj_ideal], [pt[2] for pt in traj_ideal],
    label="", linecolor=:black, linewidth=2, marker=:circle, markercolor=:red, markersize=3)

  if !isnothing(xlim) && !isnothing(ylim)
    xlims!(0, xlim)
    ylims!(0, ylim)
  end

  if !isempty(fname)
    savefig(traj_plot, fname)
  end

  return traj_plot
end

"""
    plot_trajectories(traj, traj_ideal, hc, hd, T; xlim=nothing, ylim=nothing, fname="", xlabel="\$x_1\$", ylabel="\$x_2\$", title="")

Plot the trajectories with optional axis limits and save to file if `fname` is provided.
"""
function plot_trajectories(
  traj::AbstractVector, traj_ideal::AbstractVector, hd::Real, hc::Real, T::Real;
  xlim::Union{Real,Nothing}=nothing, ylim::Union{Real,Nothing}=nothing,
  fname::String="", title::String="")
  traj_plot = plot(xlabel="\$Time [s]\$", title=title, fmt=:svg)

  traj_plot = plot!(0:hc:T, hcat(traj_ideal...)', lab=["x₁(t)" "x₂(t)"], linecolor=[:red :blue], alpha=0.6, fmt=:svg)
  traj_plot = plot!(0:hd:T, hcat(traj...)', lab=["x₁[k]" "x₂[k]"], linecolor=[:magenta :cyan], alpha=0.6, marker=:circle, markersize=2, fmt=:svg)

  if !isnothing(xlim) && !isnothing(ylim)
    xlims!(0, xlim)
    ylims!(0, ylim)
  end

  if !isempty(fname)
    savefig(traj_plot, fname)
  end

  return traj_plot
end
