import ITensorNetworks as itn
using Printf

@kwdef mutable struct TDVPProblem{State}
  state::State
  operator
  #current_time::Number = 0.0
end

state(tdvp::TDVPProblem) = tdvp.state
operator(tdvp::TDVPProblem) = tdvp.operator

function set(tdvp::TDVPProblem; state=state(tdvp), operator=operator(tdvp))
  return TDVPProblem(; state, operator)
end

function region_plan(tdvp::TDVPProblem; time_step, sweep_kwargs...) 
  return tdvp_regions(state(tdvp), time_step; sweep_kwargs...)
end

function updater(tdvp::TDVPProblem, local_tensor; region, kws...)
  local_tensor = exponentiate_updater(operator(tdvp),local_tensor; kws...)
  return tdvp, local_tensor
end

# TODO: have this callback update and print the current time
#function sweep_callback(tdvp::TDVPProblem; sweep, outputlevel, kws...)
#  if outputlevel >= 2
#    @printf("  Region %s: energy = %.12f\n",region,energy(dmrg))
#  end
#end

function applyexp(H, init_state, time_points; updater_kwargs=(;), inserter_kwargs=(;), kws...)
  init_prob = TDVPProblem(; state=copy(init_state), operator=itn.ProjTTN(H))
  time_steps = diff([0.0, time_points...])
  kwargs_array = [(; time_step=t, updater_kwargs, inserter_kwargs) for t in time_steps]
  sweep_iter = sweep_iterator(init_prob, kwargs_array)
  converged_prob = alternating_update(sweep_iter; kws...)
  return state(converged_prob)
end

