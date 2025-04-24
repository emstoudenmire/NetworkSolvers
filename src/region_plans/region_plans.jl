import Graphs: AbstractGraph, AbstractEdge, edges, dst, src, vertices
import NamedGraphs.GraphsExtensions:
  default_root_vertex, post_order_dfs_edges, post_order_dfs_vertices

function basic_forward_sweep(
  graph::AbstractGraph; nsites, root_vertex=default_root_vertex(graph), sweep_kwargs...
)
  if nsites == 1
    vertices = post_order_dfs_vertices(graph, root_vertex)
    fwd_sweep = [([v], sweep_kwargs) for v in vertices]
  elseif nsites == 2
    edges = post_order_dfs_edges(graph, root_vertex)
    fwd_sweep = [([src(e), dst(e)], sweep_kwargs) for e in edges]
  end
  return fwd_sweep
end

function basic_region_plan(args...; kws...)
  fwd_sweep = basic_forward_sweep(args...; kws...)
  # Reverse regions as well as ordering of regions
  rev_sweep = [(reverse(reg_kws[1]), reg_kws[2]) for reg_kws in reverse(fwd_sweep)]
  return [fwd_sweep..., rev_sweep...]
end
