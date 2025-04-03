import Graphs: AbstractGraph, AbstractEdge, edges, dst, src, vertices
import NamedGraphs.GraphsExtensions:
  default_root_vertex, post_order_dfs_edges, post_order_dfs_vertices

function basic_region_plan(
  graph::AbstractGraph; nsites, root_vertex=default_root_vertex(graph), sweep_kwargs...
)
  if nsites == 1
    vertices = post_order_dfs_vertices(graph, root_vertex)
    fwd_sweep = [([v], sweep_kwargs) for v in vertices]
  elseif nsites == 2
    edges = post_order_dfs_edges(graph, root_vertex)
    fwd_sweep = [([src(e), dst(e)], sweep_kwargs) for e in edges]
  end
  return [fwd_sweep..., reverse(fwd_sweep)...]
end

# TODO: add 1-site and also higher-order sweeping plans
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
  step_args = (; time_step=(-time_step / 2), dt=time_step, updater_kwargs...)

  fwd_sweep = []
  edges = post_order_dfs_edges(graph, root_vertex)
  for e in edges
    push!(fwd_sweep, ([src(e)], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))
    push!(fwd_sweep, (e, (; updater_kwargs=rev_up_args, sweep_kwargs...)))
  end
  push!(fwd_sweep, ([dst(last(edges))], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))

  # Reverse regions as well as ordering of regions:
  rev_sweep = [(reverse(rk[1]), rk[2]) for rk in reverse(fwd_sweep)]
  rev_sweep[end] = (last(rev_sweep)[1], (; updater_kwargs=step_args, sweep_kwargs...))

  return [fwd_sweep..., rev_sweep...]
end
