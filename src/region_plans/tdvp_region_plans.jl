
function tdvp_regions(graph, time_step; updater_kwargs, nsites=1, kws...)
  basic_fwd_sweep = post_order_dfs_plan(graph; nsites, kws...)

  updater_kwargs = (; nsites, time_step=(time_step/2), updater_kwargs...)

  fwd_sweep = []
  for (j, (region, region_kws)) in enumerate(basic_fwd_sweep)
    push!(fwd_sweep, (region, (; nsites, updater_kwargs, region_kws...)))
    if length(region) == 2 && j < length(basic_fwd_sweep)
      rev_kwargs = (; updater_kwargs..., time_step=(-updater_kwargs.time_step))
      push!(fwd_sweep, ([last(region)], (; updater_kwargs=rev_kwargs, region_kws...)))
    end
  end

  # Reverse regions as well as ordering of regions
  rev_sweep = [(reverse(reg_kws[1]), reg_kws[2]) for reg_kws in reverse(fwd_sweep)]
  sweep_plan = [fwd_sweep..., rev_sweep...]

  return sweep_plan
end
