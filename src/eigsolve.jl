using ITensorNetworks: ITensorNetworks
using Printf

@kwdef mutable struct EigsolveProblem{State,Operator}
  state::State
  operator::Operator
  eigenvalue::Number = Inf
end

eigenvalue(E::EigsolveProblem) = E.eigenvalue
state(E::EigsolveProblem) = E.state
operator(E::EigsolveProblem) = E.operator

function set(E::EigsolveProblem; state=state(E), operator=operator(E), eigenvalue=eigenvalue(E))
  return EigsolveProblem(; state, operator, eigenvalue)
end

function updater(E::EigsolveProblem, local_tensor; kws...)
  eigenvalue, local_tensor = eigsolve_updater(operator(E),local_tensor; kws...)
  return set(E; eigenvalue), local_tensor
end

function region_callback(E::EigsolveProblem; region, outputlevel, kws...)
  if outputlevel >= 2
    @printf("  Region %s: energy = %.12f\n",region,eigenvalue(E))
  end
  return E
end

function eigsolve(H, init_state; nsweeps, kws...)
  operator = ITensorNetworks.ProjTTN(H)
  problem = EigsolveProblem(; state=init_state, operator)
  sweep_iter = sweep_iterator(problem, fill((;),nsweeps))
  E = alternating_update(sweep_iter; kws...)
  return eigenvalue(E), state(E)
end

dmrg(args...; kws...) = eigsolve(args...; kws...)
