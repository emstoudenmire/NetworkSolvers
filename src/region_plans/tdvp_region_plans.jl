import Graphs: AbstractGraph, AbstractEdge, edges, dst, src, vertices
import NamedGraphs.GraphsExtensions:
  default_root_vertex, post_order_dfs_edges, post_order_dfs_vertices

function tdvp_regions(
  graph::AbstractGraph,
  time_step;
  root_vertex=default_root_vertex(graph),
  nsites=1,
  updater_kwargs,
  sweep_kwargs...,
)
  @assert nsites == 1
  fwd_up_args = (; time_step=(time_step / 2), dt=0.0, updater_kwargs...)
  rev_up_args = (; time_step=(-time_step / 2), dt=0.0, updater_kwargs...)
  step_args = (; time_step=(time_step / 2), dt=time_step, updater_kwargs...)

  fwd_sweep = []
  edges = post_order_dfs_edges(graph, root_vertex)
  for e in edges
    push!(fwd_sweep, ([src(e)], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))
    push!(fwd_sweep, (e, (; updater_kwargs=rev_up_args, sweep_kwargs...)))
  end
  push!(fwd_sweep, ([dst(last(edges))], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))

  # Reverse regions as well as ordering of regions:
  rev_sweep = [(reverse(rk[1]), rk[2]) for rk in reverse(fwd_sweep)]
  last_region = (first(last(rev_sweep)), (; updater_kwargs=step_args, sweep_kwargs...))
  rev_sweep = [rev_sweep[1:(end - 1)]; last_region]

  return [fwd_sweep..., rev_sweep...]
end
