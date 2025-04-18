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

function updater!(
  E::EigsolveProblem, local_tensor, region; outputlevel, solver=eigsolve_solver, kws...
)
  E.eigenvalue, local_tensor = solver(operator(E), local_tensor; kws...)
  if outputlevel >= 2
    @printf("  Region %s: energy = %.12f\n", region, eigenvalue(E))
  end
  return local_tensor
end

function eigsolve_sweep_printer(problem; outputlevel, sweep, nsweeps, kws...)
  if outputlevel >= 1
    psi = state(problem)
    print("After sweep $sweep/$nsweeps: ")
    @printf("eigenvalue=%.12f ", eigenvalue(problem))
    @printf("maxlinkdim=%d", itn.maxlinkdim(psi))
    println()
    flush(stdout)
  end
end

function eigsolve(
  H,
  init_state;
  nsweeps,
  nsites=1,
  outputlevel=0,
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  subspace_kwargs=(;),
  sweep_printer=eigsolve_sweep_printer,
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
    subspace_kwargs,
  )
  prob = alternating_update(sweep_iter; outputlevel, sweep_printer, kws...)
  return eigenvalue(prob), state(prob)
end

dmrg(args...; kws...) = eigsolve(args...; kws...)
