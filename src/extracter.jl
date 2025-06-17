import ConstructionBase: setproperties

function extracter(problem, region_iterator; kws...)
  region = current_region(region_iterator)
  psi = itn.orthogonalize(state(problem), region)
  local_state = prod(psi[v] for v in region)
  return setproperties(problem; state=psi), local_state
end
