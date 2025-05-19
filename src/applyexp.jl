import ITensorNetworks as itn
using Printf

@kwdef mutable struct ApplyExpProblem{State}
  state::State
  operator
  current_time::Number = 0.0
end

state(tdvp::ApplyExpProblem) = tdvp.state
operator(tdvp::ApplyExpProblem) = tdvp.operator
current_time(tdvp::ApplyExpProblem) = tdvp.current_time

function set!(
  T::ApplyExpProblem; state=state(T), operator=operator(T), current_time=current_time(T)
)
  T.state = state
  T.operator = operator
  T.current_time = current_time
end

function region_plan(tdvp::ApplyExpProblem; nsites, time_step, sweep_kwargs...)
  return tdvp_regions(state(tdvp), time_step; nsites, sweep_kwargs...)
end

function updater!(
  T::ApplyExpProblem,
  local_state,
  region_iterator;
  nsites,
  time_step,
  solver=runge_kutta_solver,
  outputlevel,
  kws...,
)
  local_state, info = solver(x->optimal_map(operator(T), x), time_step, local_state; kws...)

  if nsites==1
    curr_reg = current_region(region_iterator)
    next_reg = next_region(region_iterator)
    if !isnothing(next_reg) && next_reg != curr_reg
      seq = itn.edge_sequence_between_regions(state(T), curr_reg, next_reg)
      next_edge = first(seq)
      v1, v2 = itn.src(next_edge), itn.dst(next_edge)
      psi = copy(state(T))
      psi[v1], R = qr(local_state, uniqueinds(local_state, psi[v2]))
      shifted_operator = itn.position(operator(T), psi, itn.NamedEdge(v1=>v2))
      R_t, _ = solver(x->optimal_map(shifted_operator, x), -time_step, R; kws...)
      local_state = psi[v1]*R_t
    end
  end

  if is_last_region(region_iterator)
    # TODO: move this to a new "applyexp_sweep_printer" function
    T.current_time += 2*abs(time_step)  # currently assuming second-order method
    if outputlevel >= 1
      @printf("  Current time = %s, ", current_time(T))
      @printf("maxlinkdim=%d", itn.maxlinkdim(state(T)))
      println()
      flush(stdout)
    end
  end
  return local_state
end

function applyexp(
  init_prob,
  exponents;
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  outputlevel=0,
  nsites=1,
  subspace_kwargs=(;),
  tdvp_order=4,
  kws...,
)
  time_steps = diff([0.0, exponents...])[2:end]
  sweep_kws = (;
    outputlevel,
    extracter_kwargs,
    inserter_kwargs,
    nsites,
    subspace_kwargs,
    tdvp_order,
    updater_kwargs,
  )
  kws_array = [(; sweep_kws..., time_step=t) for t in time_steps]
  sweep_iter = sweep_iterator(init_prob, kws_array)
  converged_prob = sweep_solve(sweep_iter; outputlevel, kws...)
  return state(converged_prob)
end

function applyexp(H, init_state, exponents; kws...)
  init_prob = ApplyExpProblem(;
    state=permute_indices(init_state), operator=itn.ProjTTN(permute_indices(H))
  )
  return applyexp(init_prob, exponents; kws...)
end

function tdvp(H, init_state, time_points; time_angle=0.0, kws...)
  z = exp(-im*time_angle)
  exponents = [(-im*z)*t for t in time_points]
  return applyexp(H, init_state, exponents; kws...)
end
