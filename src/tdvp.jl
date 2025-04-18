import ITensorNetworks as itn
using Printf

@kwdef mutable struct TDVPProblem{State}
  state::State
  operator
  current_time::Number = 0.0
end

state(tdvp::TDVPProblem) = tdvp.state
operator(tdvp::TDVPProblem) = tdvp.operator
current_time(tdvp::TDVPProblem) = tdvp.current_time

function region_plan(tdvp::TDVPProblem; nsites, time_step, sweep_kwargs...)
  return tdvp_regions(state(tdvp), time_step; nsites, sweep_kwargs...)
end

function updater!(
  T::TDVPProblem,
  local_tensor,
  region;
  time_step,
  dt,
  solver=exponentiate_solver,
  outputlevel,
  kws...,
)
  local_tensor, info = solver(operator(T), time_step, local_tensor; kws...)
  T.current_time += dt
  if outputlevel >= 2 && abs(dt) > 0.0
    @printf("  Current time = %s\n", current_time(T))
  end
  return local_tensor
end

function applyexp(
  H,
  init_state,
  time_points;
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  outputlevel=0,
  nsites=1,
  subspace_kwargs=(;),
  kws...,
)
  init_prob = TDVPProblem(; state=copy(init_state), operator=itn.ProjTTN(H))
  time_steps = diff([0.0, time_points...])[2:end]
  sweep_kws = (;
    outputlevel, extracter_kwargs, inserter_kwargs, nsites, subspace_kwargs, updater_kwargs
  )
  kws_array = [(; sweep_kws..., time_step=t) for t in time_steps]
  sweep_iter = sweep_iterator(init_prob, kws_array)
  converged_prob = alternating_update(sweep_iter; outputlevel, kws...)
  return state(converged_prob)
end
