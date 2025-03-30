import ITensors as it
import ITensorNetworks as itn
using NamedGraphs: NamedGraphs
using Graphs: Graphs

default_maxdim() = typemax(Int)
default_mindim() = 1
default_cutoff() = 0.0

function extracter(problem; region, kws...)
  psi = copy(state(problem))
  if isa(region, Graphs.AbstractEdge)
    vsrc, vdst = Graphs.src(region), Graphs.dst(region)
    psi = itn.orthogonalize(psi, vsrc)
    left_inds = it.uniqueinds(psi[vsrc], psi[vdst])
    lefttags = it.tags(psi, region)
    righttags = it.tags(psi, region)
    U, S, V = it.svd(psi[vsrc], left_inds; lefttags, righttags)
    psi[vsrc] = U
    local_tensor = S * V
  else
    psi = itn.orthogonalize(psi, region)
    local_tensor = prod(psi[v] for v in region)
  end
  op = itn.position(operator(problem), psi, region)
  return set(problem; state=psi, operator=op), local_tensor
end

function inserter(
  problem,
  local_tensor,
  region;
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  cutoff=default_cutoff(),
  normalize=false,
  kws...,
)
  psi = copy(state(problem))
  v = last(region)
  if length(region) == 2
    e = NamedGraphs.edgetype(psi)(first(region), last(region))
    indsTe = it.inds(psi[first(region)])
    tags = it.tags(psi, e)
    U, C, _ = it.factorize(local_tensor, indsTe; tags, maxdim, mindim, cutoff)
    psi[first(region)] = U
  elseif length(region) == 1
    C = local_tensor
  else
    error("Only length==2 or length==1 regions currently supported")
  end
  psi[v] = C
  psi = itn.set_ortho_region(psi, [v])
  normalize && (psi[v] /= norm(psi[v]))
  return set(problem; state=psi)
end

function inserter(problem, local_tensor, region::NamedGraphs.NamedEdge; kws...)
  #TODO: potential bug / workaround in this function
  #      TreeTensorNetwork type does not allow copying
  #      if `is_tree(ttn)` is false. This is probably
  #      too restrictive, since copying is a low-level operation.
  #      Would like to write this function more like code above,
  #      without mutating the state referenced by `problem`.
  #psi = copy(state(problem))
  psi = state(problem)
  psi[Graphs.dst(region)] *= local_tensor
  psi = itn.set_ortho_region(psi, [Graphs.dst(region)])
  return set(problem; state=psi)
end
