import Graphs: AbstractGraph, edges, dst, src, vertices

function basic_path_regions(g::AbstractGraph; sweep_kwargs...)
  fwd_sweep = [([src(e), dst(e)], sweep_kwargs) for e in edges(g)]
  return [fwd_sweep..., reverse(fwd_sweep)...]
end

function tdvp_regions(g::AbstractGraph, time_step; updater_kwargs, sweep_kwargs...)
  fwd_up_args = (; time=(time_step/2), updater_kwargs...)
  rev_up_args = (; time=(-time_step/2), updater_kwargs...)

  fwd_sweep = []
  for e in edges(g)
    push!(fwd_sweep,([src(e)], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))
    push!(fwd_sweep,(e, (; updater_kwargs=rev_up_args, sweep_kwargs...)))
  end
  push!(fwd_sweep,([dst(last(edges(g)))], (; updater_kwargs=fwd_up_args, sweep_kwargs...)))

  # Reverse regions as well as ordering of regions:
  rev_sweep = [(reverse(rk[1]),rk[2]) for rk in reverse(fwd_sweep)]

  return [fwd_sweep..., rev_sweep...]
end
