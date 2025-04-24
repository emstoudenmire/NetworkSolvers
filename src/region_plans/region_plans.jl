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
