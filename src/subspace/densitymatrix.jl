using ITensors: pause
import ITensorNetworks as itn

function subspace_expand!(
  ::Backend"densitymatrix",
  problem::EigsolveProblem,
  local_tensor,
  region;
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  north_pass=2,
  expansion_factor=1.5,
  kws...,
)
  isa(region, AbstractEdge) && return local_tensor
  psi = state(problem)

  prev_vertex_set = setdiff(itn.pos(operator(problem)), region)
  (length(prev_vertex_set) != 1) && return local_tensor
  prev_vertex = only(prev_vertex_set)
  A = psi[prev_vertex]

  next_vertices = filter(v -> (it.hascommoninds(psi[v], A)), region)
  isempty(next_vertices) && return local_tensor
  next_vertex = only(next_vertices)
  C = psi[next_vertex]

  a = commonind(A, C)
  isnothing(a) && return local_tensor

  basis_size = prod(dim.(uniqueinds(A, C)))
  expand_maxdim = min(maxdim, ceil(Int, (expansion_factor-1) * dim(a)))
  expand_maxdim = min(basis_size-dim(a), expand_maxdim)
  expand_maxdim <= 0 && return local_tensor

  env = itn.environments(operator(problem))
  incident_edges = itn.incident_edges(operator(problem))
  prev_env_edges = filter(e->(src(e) ∉ region && dst(e) ∉ region), incident_edges)
  E = isempty(prev_env_edges) ? ITensor(1.0) : env[only(prev_env_edges)]

  H = itn.operator(operator(problem))
  sqrt_rho = E*A*H[prev_vertex]
  conj_proj_A(T) = (T - prime(A)*(dag(prime(A))*T))
  for pass in 1:north_pass
    sqrt_rho = conj_proj_A(sqrt_rho)
  end
  rho = sqrt_rho * dag(noprime(sqrt_rho))
  D, U = eigen(rho; cutoff, maxdim=expand_maxdim, mindim, ishermitian=true)

  ###TODO: do we need this?
  Uproj(T) = (T - prime(A, a)*(dag(prime(A, a))*T))
  U = Uproj(U)
  ###

  if norm(dag(U)*A) > 1E-10
    @printf("Warning: |U*A| = %.3E in subspace expansion\n", norm(dag(U)*A))
    return local_tensor
  end

  Ax, ax = directsum(A=>a, U=>commonind(U, D))

  expander = dag(Ax) * A
  psi[prev_vertex] = Ax
  psi[next_vertex] = expander * C

  # TODO: avoid computing local tensor twice
  #       while also handling AbstractEdge region case
  local_tensor = prod(psi[v] for v in region)

  return local_tensor
end
