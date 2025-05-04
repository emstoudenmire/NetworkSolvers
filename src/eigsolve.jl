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

function set!(
  E::EigsolveProblem; state=E.state, operator=E.operator, eigenvalue=E.eigenvalue
)
  E = EigsolveProblem(; state, operator, eigenvalue)
end

function updater!(
  E::EigsolveProblem,
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

function eigsolve_sweep_printer(problem; outputlevel, sweep, nsweeps, kws...)
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
    flush(stdout)
  end
end

function eigsolve(
  init_prob;
  nsweeps,
  nsites=1,
  outputlevel=0,
  extracter_kwargs=(;),
  updater_kwargs=(;),
  truncation_kwargs=(;),
  subspace_kwargs=(; algorithm="densitymatrix"),
  sweep_printer=eigsolve_sweep_printer,
  kws...,
)
  sweep_iter = sweep_iterator(
    init_prob,
    nsweeps;
    nsites,
    outputlevel,
    extracter_kwargs,
    updater_kwargs,
    truncation_kwargs,
    subspace_kwargs,
  )
  prob = alternating_update(sweep_iter; outputlevel, sweep_printer, kws...)
  return eigenvalue(prob), state(prob)
end

function eigsolve(H, init_state; kws...)
  init_prob = EigsolveProblem(;
    state=permute_indices(init_state), operator=itn.ProjTTN(permute_indices(H))
  )
  return eigsolve(init_prob; kws...)
end

dmrg(args...; kws...) = eigsolve(args...; kws...)
