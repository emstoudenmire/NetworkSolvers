import Graphs: AbstractGraph, edges, dst, src, vertices

function basic_path_regions(g::AbstractGraph; sweep_kwargs...)
  fwd_sweep = [([src(e), dst(e)], sweep_kwargs) for e in edges(g)]
  return [fwd_sweep..., reverse(fwd_sweep)...]
end

# TODO: does this region plan show that the idea of "sweep kwargs"
# is misleading? Because this makes two sweeps, and we can have four sweeps etc.
# Maybe "iteration kwargs" ?
function tdvp_regions(g::AbstractGraph, time_step; sweep_kwargs...)
  verts = vertices(g)
  fwd_sweep = []
  for e in edges(g)
    push!(fwd_sweep,([src(e)], (; time_step=(time_step / 2), sweep_kwargs...)))
    push!(fwd_sweep,(e, (; time_step=(-time_step / 2), sweep_kwargs...)))
  end
  push!(fwd_sweep,([dst(last(edges(g)))], (; time_step=(time_step / 2), sweep_kwargs...)))
  # Reverse regions as well as ordering of regions:
  rev_sweep = [(reverse(rk[1]),rk[2]) for rk in reverse(fwd_sweep)]
  return [fwd_sweep..., rev_sweep...]
end
