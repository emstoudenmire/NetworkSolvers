
function extracter!(problem, region_iterator; kws...)
  region = current_region(region_iterator)
  psi = itn.orthogonalize(state(problem), region)
  set!(problem; state=psi)
  local_state = prod(psi[v] for v in region)
  return local_state
end
