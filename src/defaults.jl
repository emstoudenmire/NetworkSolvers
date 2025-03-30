import ITensors as it
import ITensorNetworks as itn
using NamedGraphs: NamedGraphs
using Graphs: Graphs

default_maxdim() = typemax(Int)
default_mindim() = 1
default_cutoff() = 0.0

function extracter!(problem, region; kws...)
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
  problem.operator = itn.position(operator(problem), psi, region)
  problem.state = psi
  return local_tensor
end

function inserter!(
  problem,
  local_tensor,
  region;
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  cutoff=default_cutoff(),
  normalize=false,
  kws...,
)
  psi = state(problem)
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
  problem.state = psi
  return
end

function inserter!(problem, local_tensor, region::NamedGraphs.NamedEdge; kws...)
  psi = state(problem)
  psi[Graphs.dst(region)] *= local_tensor
  problem.state = itn.set_ortho_region(psi, [Graphs.dst(region)])
  return
end

region_printer(problem, local_tensor, region; kws...) = nothing
