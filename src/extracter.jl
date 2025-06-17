import ConstructionBase: setproperties

function extracter(problem, region_iterator; sweep, kws...)
  region = current_region(region_iterator)
  psi = itn.orthogonalize(state(problem), region)
  local_state = prod(psi[v] for v in region)
  problem, local_state = subspace_expand(
    problem, local_state, region_iterator; sweep, kws...
  )
  shifted_operator = itn.position(operator(problem), psi, region)
  return setproperties(problem; state=psi, operator=shifted_operator), local_state
end
