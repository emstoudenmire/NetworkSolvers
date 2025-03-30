import ITensorNetworks as itn
using Printf

@kwdef mutable struct EigsolveProblem{State,Operator}
  state::State
  operator::Operator
  eigenvalue::Number = Inf
end

eigenvalue(E::EigsolveProblem) = E.eigenvalue
state(E::EigsolveProblem) = E.state
operator(E::EigsolveProblem) = E.operator

function set(
  E::EigsolveProblem; state=state(E), operator=operator(E), eigenvalue=eigenvalue(E)
)
  return EigsolveProblem(; state, operator, eigenvalue)
end

function updater(E::EigsolveProblem, local_tensor; region, kws...)
  eigenvalue, local_tensor = eigsolve_updater(operator(E), local_tensor; kws...)
  return set(E; eigenvalue), local_tensor
end

function region_callback(E::EigsolveProblem; region, outputlevel, kws...)
  if outputlevel >= 2
    @printf("  Region %s: energy = %.12f\n", region, eigenvalue(E))
  end
  return E
end

function eigsolve(
  H, init_state; nsweeps, nsites=2, updater_kwargs=(;), inserter_kwargs=(;), kws...
)
  init_prob = EigsolveProblem(; state=copy(init_state), operator=itn.ProjTTN(H))
  kwargs_array = [(; nsites, sweep=sw, updater_kwargs, inserter_kwargs) for sw in 1:nsweeps]
  sweep_iter = sweep_iterator(init_prob, kwargs_array)
  converged_prob = alternating_update(sweep_iter; kws...)
  return eigenvalue(converged_prob), state(converged_prob)
end

dmrg(args...; kws...) = eigsolve(args...; kws...)
