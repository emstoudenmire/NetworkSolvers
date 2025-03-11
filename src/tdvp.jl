using ITensorNetworks: ITensorNetworks
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
  return tdvp_regions(state(tdvp),time_step; sweep_kwargs...)
end

function updater(tdvp::TDVPProblem, local_tensor; exponentiate_kwargs, time_step, kws...)
  local_tensor = exponentiate_updater(operator(tdvp),local_tensor; time_step, exponentiate_kwargs...)
  return tdvp, local_tensor
end

# TODO: have this callback update and print the current time
#function sweep_callback(tdvp::TDVPProblem; sweep, outputlevel, kws...)
#  if outputlevel >= 2
#    @printf("  Region %s: energy = %.12f\n",region,energy(dmrg))
#  end
#end

"""
  tdvp(H, init_state; time_steps, nsites=1, kwargs...)

Accepts an array or iterable of time steps.
"""
function tdvp(H, init_state; time_steps, nsites=1, exponentiate_kwargs=(;), inserter_kwargs=(;), kws...)
  operator = ITensorNetworks.ProjTTN(H)
  kwargs_array = [(; time_step=t, exponentiate_kwargs, inserter_kwargs) for t in time_steps]
  problem = TDVPProblem(;state=init_state,operator)
  sweep_iter = sweep_iterator(problem, kwargs_array)
  tdvp_prob = alternating_update(sweep_iter; kws...)
  return state(tdvp_prob)
end

"""
  tdvp(H, init_state, times; kws...)

Accepts an array of time points.
"""
function tdvp(H, init_state, times; kws...)
  time_steps = Vector{eltype(times)}(undef,length(times))
  time_steps[1] = times[1]
  for j in 2:length(times)
    time_steps[j] = times[j] - times[j - 1]
  end
  return tdvp(H, init_state; time_steps, kws...)
end

"""
  tdvp(H, init_state, times; kws...)

Accepts a total time and time step.
"""
function tdvp(H, init_state, total_time::Number, time_step::Number; kws...)
  nsteps = convert(Int, ceil(abs(total_time / time_step)))
  @assert (nsteps * time_step ≈ total_time)
  return tdvp(H, init_state; time_steps=fill(time_step,nsteps), kws...)
end
