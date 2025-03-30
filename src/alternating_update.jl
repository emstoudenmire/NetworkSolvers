
function region_callback(problem; kws...)
  return problem
end

function sweep_callback(problem; outputlevel, sweep, kws...)
  if outputlevel >= 1
    println("Done with sweep $sweep")
    #print(" cpu_time=", round(sweep_time; digits=3))
    #println()
    flush(stdout)
  end
  return problem
end

function alternating_update(
  sweep_iterator;
  outputlevel=0,
  kwargs...,
)
  for (sweep, region_iter) in enumerate(sweep_iterator)
    prob = problem(region_iter)
    for (region, region_kwargs) in region_tuples(region_iter)
      # TODO: maybe just don't let region callback change the problem!
      prob = region_callback(prob;
        region,
        nsweeps=length(sweep_iterator),
        outputlevel,
        sweep,
        region_kwargs...,
        kwargs...,
      )
      #set(sweep_iterator; problem) #TODO
    end
    prob = sweep_callback(prob; nsweeps=length(sweep_iterator), outputlevel, sweep, kwargs...)
    #TODO need to update `problem` reference held within region iterator(s)
  end
  return problem(last(sweep_iterator))
end

