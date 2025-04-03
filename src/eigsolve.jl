import ITensorNetworks as itn
import ITensors as it
using Printf

@kwdef mutable struct EigsolveProblem{State,Operator}
  state::State
  operator::Operator
  eigenvalue::Number = Inf
end

eigenvalue(E::EigsolveProblem) = E.eigenvalue
state(E::EigsolveProblem) = E.state
operator(E::EigsolveProblem) = E.operator

function updater!(E::EigsolveProblem, local_tensor, region; outputlevel, kws...)
  E.eigenvalue, local_tensor = eigsolve_updater(operator(E), local_tensor; kws...)
  if outputlevel >= 2
    @printf("  Region %s: energy = %.12f\n", region, eigenvalue(E))
  end
  return local_tensor
end

function eigsolve(
  H,
  init_state;
  nsweeps,
  nsites=2,
  outputlevel=0,
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  kws...,
)
  init_prob = EigsolveProblem(; state=copy(init_state), operator=itn.ProjTTN(H))
  sweep_iter = sweep_iterator(
    init_prob,
    nsweeps;
    nsites,
    outputlevel,
    extracter_kwargs,
    updater_kwargs,
    inserter_kwargs,
  )
  prob = alternating_update(sweep_iter; outputlevel, kws...)
  return eigenvalue(prob), state(prob)
end

dmrg(args...; kws...) = eigsolve(args...; kws...)
