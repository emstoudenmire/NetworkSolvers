
function inserter!(
  problem,
  local_tensor,
  region_iterator;
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  normalize=false,
  set_orthogonal_region=true,
  sweep,
  kws...,
)
  cutoff = get_or_last(cutoff, sweep)
  mindim = get_or_last(mindim, sweep)
  maxdim = get_or_last(maxdim, sweep)

  region = current_region(region_iterator)
  psi = state(problem)
  if region isa ng.NamedEdge
    psi = state(problem)
    #Is there a world in which this modifies the graph structure, if not setindex also should be used here for efficiency
    psi[Graphs.dst(region)] *= local_tensor
    psi = itn.set_ortho_region(psi, [Graphs.dst(region)])
    set!(problem; state=psi)
    return nothing
  elseif length(region) == 2
    e = ng.edgetype(psi)(first(region), last(region))
    indsTe = it.inds(psi[first(region)])
    tags = it.tags(psi, e)
    U, C, _ = it.factorize(local_tensor, indsTe; tags, maxdim, mindim, cutoff)
    itn.setindex_preserve_graph!(psi, U, first(region))
  elseif length(region) == 1
    C = local_tensor
  else
    error("Region of length $(length(region)) not currently supported")
  end
  v = last(region)
  itn.setindex_preserve_graph!(psi, C, v)
  psi = set_orthogonal_region ? itn.set_ortho_region(psi, [v]) : psi
  normalize && itn.setindex_preserve_graph!(psi, psi[v] / norm(psi[v]), v)
  set!(problem; state=psi)
  return nothing
end
