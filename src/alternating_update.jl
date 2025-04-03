
default_region_callback(problem; kws...) = nothing

default_sweep_callback(problem; kws...) = nothing

function default_sweep_printer(problem; outputlevel, sweep, kws...)
  if outputlevel >= 1
    println("Done with sweep $sweep")
    #println(" CPU time=", round(sweep_time; digits=3))
    flush(stdout)
  end
end

function alternating_update(
  sweep_iterator;
  outputlevel=0,
  region_callback=default_region_callback,
  sweep_callback=default_sweep_callback,
  sweep_printer=default_sweep_printer,
  kwargs...,
)
  for (sweep, region_iter) in enumerate(sweep_iterator)
    prob = problem(region_iter)
    for (region, region_kwargs) in region_tuples(region_iter)
      region_callback(
        prob;
        nsweeps=length(sweep_iterator),
        outputlevel,
        region,
        region_kwargs,
        sweep,
        kwargs...,
      )
    end
    sweep_callback(prob; nsweeps=length(sweep_iterator), outputlevel, sweep, kwargs...)
    sweep_printer(prob; nsweeps=length(sweep_iterator), outputlevel, sweep, kwargs...)
  end
  return problem(last(sweep_iterator))
end
