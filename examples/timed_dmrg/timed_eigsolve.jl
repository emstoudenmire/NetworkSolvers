import ITensorNetworks as itn
using Printf
import NetworkSolvers: current_region, region_iterator_action!, extracter!, updater!, inserter!, prepare_subspace!, eigenvalue, state, operator, eigsolve, eigsolve_solver, permute_indices

@kwdef mutable struct ProblemTimings
  extracter_time::Float64 = 0.0
  subspace_time::Float64 = 0.0
  updater_time::Float64 = 0.0
  inserter_time::Float64 = 0.0
end

function print(T::ProblemTimings)
  @printf("  Extracter time = %.3f s\n", T.extracter_time)
  @printf("  Subspace time = %.3f s\n", T.subspace_time)
  @printf("  Updater time = %.3f s\n", T.updater_time)
  @printf("  Inserter time = %.3f s\n", T.inserter_time)
end

@kwdef mutable struct TimedEigsolveProblem{State,Operator}
  state::State
  operator::Operator
  eigenvalue::Number = Inf
  timings::ProblemTimings = ProblemTimings()
end

eigenvalue(E::TimedEigsolveProblem) = E.eigenvalue
state(E::TimedEigsolveProblem) = E.state
operator(E::TimedEigsolveProblem) = E.operator

function updater!(
  E::TimedEigsolveProblem,
  local_tensor,
  region_iterator;
  outputlevel,
  solver=eigsolve_solver,
  kws...,
)
  E.eigenvalue, local_tensor = solver(operator(E), local_tensor; kws...)
  if outputlevel >= 2
    @printf("  Region %s: energy = %.12f\n", current_region(region_iterator), eigenvalue(E))
  end
  return local_tensor
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
  problem.timings.extracter_time += @elapsed begin
    local_tensor = extracter!(problem, region_iterator; extracter_kwargs..., kwargs...)
  end
  problem.timings.subspace_time += @elapsed begin
    local_tensor = prepare_subspace!(
      problem, local_tensor, region_iterator; subspace_kwargs..., sweep, kwargs...
    )
  end
  problem.timings.updater_time += @elapsed begin
    local_tensor = updater!(problem, local_tensor, region_iterator; updater_kwargs..., kwargs...)
  end
  problem.timings.inserter_time += @elapsed begin
    inserter!(problem, local_tensor, region_iterator; sweep, truncation_kwargs..., kwargs...)
  end
  return nothing
end

function timed_eigsolve_sweep_printer(problem; outputlevel, sweep, nsweeps, kws...)
  if outputlevel >= 1
    psi = state(problem)
    if nsweeps >= 10
      @printf("After sweep %02d/%d ", sweep, nsweeps)
    else
      @printf("After sweep %d/%d ", sweep, nsweeps)
    end
    @printf("eigenvalue=%.12f ", eigenvalue(problem))
    @printf("maxlinkdim=%d", itn.maxlinkdim(psi))
    println()
    print(problem.timings)
    problem.timings = ProblemTimings()
    println()
    flush(stdout)
  end
end

function timed_eigsolve(
  H,
  init_state;
  sweep_printer=timed_eigsolve_sweep_printer,
  kws...,
)
  init_prob = TimedEigsolveProblem(;
    state=permute_indices(init_state), operator=itn.ProjTTN(permute_indices(H))
  )
  return eigsolve(init_prob; sweep_printer, kws...)
end

timed_dmrg(args...; kws...) = timed_eigsolve(args...; kws...)
