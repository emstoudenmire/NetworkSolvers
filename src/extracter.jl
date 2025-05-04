
function extracter!(problem, region_iterator; kws...)
  region = current_region(region_iterator)
  if isa(region, Graphs.AbstractEdge)
    vsrc, vdst = Graphs.src(region), Graphs.dst(region)
    psi = itn.orthogonalize(state(problem), vsrc)
    left_inds = it.uniqueinds(psi[vsrc], psi[vdst])
    lefttags = it.tags(psi, region)
    righttags = it.tags(psi, region)
    U, S, V = it.svd(psi[vsrc], left_inds; lefttags, righttags)
    psi[vsrc] = U
    local_tensor = S * V
  else
    psi = itn.orthogonalize(state(problem), region)
    local_tensor = prod(psi[v] for v in region)
  end
  set!(problem; state=psi)
  return local_tensor
end
