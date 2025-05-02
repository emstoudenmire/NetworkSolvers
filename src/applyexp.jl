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
  region_iterator;
  time_step,
  solver=exponentiate_solver,
  outputlevel,
  kws...,
)
  local_tensor, info = solver(operator(T), time_step, local_tensor; kws...)
  if outputlevel >= 2 && is_last_region(region_iterator)
    T.current_time += 2*time_step  # currently assuming second-order method
    @printf("  Current time = %s\n", current_time(T))
  end
  return local_tensor
end

function applyexp(
  H,
  init_state,
  exponents;
  extracter_kwargs=(;),
  updater_kwargs=(;),
  truncation_kwargs=(;),
  outputlevel=0,
  nsites=1,
  subspace_kwargs=(;),
  kws...,
)
  H = permute_indices(H)
  init_state = permute_indices(init_state)
  init_prob = TDVPProblem(; state=copy(init_state), operator=itn.ProjTTN(H))
  time_steps = diff([0.0, exponents...])[2:end]
  sweep_kws = (;
    outputlevel,
    extracter_kwargs,
    truncation_kwargs,
    nsites,
    subspace_kwargs,
    updater_kwargs,
  )
  kws_array = [(; sweep_kws..., time_step=t) for t in time_steps]
  sweep_iter = sweep_iterator(init_prob, kws_array)
  converged_prob = alternating_update(sweep_iter; outputlevel, kws...)
  return state(converged_prob)
end

function tdvp(H, init_state, time_points; time_angle=0.0, kws...)
  z = exp(-im*time_angle)
  exponents = [(-im*z)*t for t in time_points]
  return applyexp(H, init_state, exponents; kws...)
end
