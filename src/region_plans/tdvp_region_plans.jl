import Graphs: AbstractGraph, AbstractEdge, edges, dst, src, vertices
import NamedGraphs.GraphsExtensions:
  default_root_vertex, post_order_dfs_edges, post_order_dfs_vertices

function tdvp_regions(graph::AbstractGraph, time_step; updater_kwargs, kws...)
  basic_fwd_sweep = basic_forward_sweep(graph; kws...)

  fwd_update = (; time_step=(+time_step/2), updater_kwargs...)
  rev_update = (; time_step=(-time_step/2), updater_kwargs...)

  fwd_sweep = []
  for (j, (region, region_kws)) in enumerate(basic_fwd_sweep)
    push!(fwd_sweep, (region, (; updater_kwargs=fwd_update, region_kws...)))
    # Put in reverse step except at end of forward sweep
    if j < length(basic_fwd_sweep)
      if length(region) == 1
        next_region = first(basic_fwd_sweep[j + 1])
        rev_region = NamedEdge(only(region), only(next_region))
      elseif length(region) == 2
        rev_region = [last(region)]
      else
        error("TDVP currently does not support regions of length = $(length(region))")
      end
      push!(fwd_sweep, (rev_region, (; updater_kwargs=rev_update, region_kws...)))
    end
  end

  # Reverse regions as well as ordering of regions
  rev_sweep = [(reverse(reg_kws[1]), reg_kws[2]) for reg_kws in reverse(fwd_sweep)]
  sweep_plan = [fwd_sweep..., rev_sweep...]

  return sweep_plan
end
