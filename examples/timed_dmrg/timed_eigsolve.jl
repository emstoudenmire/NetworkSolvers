import ITensorNetworks as itn
using Printf
import NetworkSolvers:
  EigsolveProblem,
  current_region,
  region_iterator_action!,
  extracter!,
  updater!,
  inserter!,
  prepare_subspace!,
  eigenvalue,
  set!,
  state,
  operator,
  eigsolve,
  eigsolve_solver,
  eigsolve_sweep_printer,
  permute_indices

@kwdef mutable struct TimedEigsolveProblem{State,Operator}
  eigprob::EigsolveProblem{State,Operator}
  extracter_time::Float64 = 0.0
  subspace_time::Float64 = 0.0
  updater_time::Float64 = 0.0
  inserter_time::Float64 = 0.0
end

reset_timings(T::TimedEigsolveProblem) = TimedEigsolveProblem(; eigprob=T.eigprob)
eigenvalue(E::TimedEigsolveProblem) = eigenvalue(E.eigprob)
state(E::TimedEigsolveProblem) = state(E.eigprob)
operator(E::TimedEigsolveProblem) = operator(E.eigprob)

function set!(
  T::TimedEigsolveProblem; state=state(T), operator=operator(T), eigenvalue=eigenvalue(T)
)
  T.eigprob = EigsolveProblem(; state, operator, eigenvalue)
end

function updater!(E::TimedEigsolveProblem, local_tensor, region_iterator; kws...)
  return updater!(E.eigprob, local_tensor, region_iterator; kws...)
end

function region_iterator_action!(
  problem::TimedEigsolveProblem,
  region_iterator;
  extracter_kwargs=(;),
  subspace_kwargs=(;),
  updater_kwargs=(;),
  truncation_kwargs=(;),
  sweep,
  kwargs...,
)
  problem.extracter_time += @elapsed begin
    local_tensor = extracter!(problem, region_iterator; extracter_kwargs..., kwargs...)
  end
  problem.subspace_time += @elapsed begin
    local_tensor = prepare_subspace!(
      problem, local_tensor, region_iterator; subspace_kwargs..., sweep, kwargs...
    )
  end
  problem.updater_time += @elapsed begin
    local_tensor = updater!(
      problem, local_tensor, region_iterator; updater_kwargs..., kwargs...
    )
  end
  problem.inserter_time += @elapsed begin
    inserter!(
      problem, local_tensor, region_iterator; sweep, truncation_kwargs..., kwargs...
    )
  end
  return nothing
end

function timed_eigsolve_sweep_printer(problem; outputlevel, kws...)
  eigsolve_sweep_printer(problem.eigprob; outputlevel, kws...)
  if outputlevel >= 1
    @printf("  Extracter time = %.3f s\n", problem.extracter_time)
    @printf("  Subspace time = %.3f s\n", problem.subspace_time)
    @printf("  Updater time = %.3f s\n", problem.updater_time)
    @printf("  Inserter time = %.3f s\n", problem.inserter_time)
    problem = reset_timings(problem)
    println()
    flush(stdout)
  end
end

function timed_eigsolve(H, ψ0; sweep_printer=timed_eigsolve_sweep_printer, kws...)
  eigprob = EigsolveProblem(;
    state=permute_indices(ψ0), operator=itn.ProjTTN(permute_indices(H))
  )
  return eigsolve(TimedEigsolveProblem(; eigprob); sweep_printer, kws...)
end

timed_dmrg(args...; kws...) = timed_eigsolve(args...; kws...)
