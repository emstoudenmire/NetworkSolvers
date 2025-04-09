import ITensorNetworks as itn

function subspace_expand!(
  ::Backend"densitymatrix",
  problem::EigsolveProblem,
  local_tensor,
  region;
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  north_pass=1,
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

  envs = itn.environments(operator(problem))
  H = itn.operator(operator(problem))
  sqrt_rho = A
  for e in itn.incident_edges(operator(problem))
    (src(e) ∈ region || dst(e) ∈ region) && continue
    sqrt_rho *= envs[e]
  end
  sqrt_rho *= H[prev_vertex]

  conj_proj_A(T) = (T - prime(A)*(dag(prime(A))*T))
  for pass in 1:north_pass
    sqrt_rho = conj_proj_A(sqrt_rho)
  end
  rho = sqrt_rho * dag(noprime(sqrt_rho))
  D, U = eigen(rho; cutoff, maxdim=expand_maxdim, mindim, ishermitian=true)

  Uproj(T) = (T - prime(A, a)*(dag(prime(A, a))*T))
  for pass in 1:north_pass
    U = Uproj(U)
  end
  if norm(dag(U)*A) > 1E-10
    @printf("Warning: |U*A| = %.3E in subspace expansion\n", norm(dag(U)*A))
    return local_tensor
  end

  Ax, ax = directsum(A=>a, U=>commonind(U, D))
  expander = dag(Ax) * A
  psi[prev_vertex] = Ax
  psi[next_vertex] = expander * C
  local_tensor = expander*local_tensor

  return local_tensor
end
