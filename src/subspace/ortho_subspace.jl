using ITensors

expand_space(χ::Integer, expansion_factor) = max(χ + 1, floor(Int, expansion_factor * χ))

function expand_space(χs::Vector{<:Pair}, expansion_factor)
  return [q => expand_space(d, expansion_factor) for (q, d) in χs]
end

#
# Alternative idea for "ortho" method:
#  - Just make a random tensor with `basis_inds` on both sides
#    (Kind of like a random density matrix.)
#  - Symmetrize to make it Hermitian PSD.
#  - Then do eigenvalue decomposition to get U at desired size.
#  - Finally, project out space of A from U.
#

function subspace_expand!(
  ::Backend"ortho",
  problem::EigsolveProblem,
  local_tensor,
  region;
  prev_region,
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  expansion_factor=1.2,
  kws...,
)
  if isnothing(prev_region) || isa(region, AbstractEdge)
    return local_tensor
  end

  prev_vertex_set = setdiff(prev_region, region)
  (length(prev_vertex_set) != 1) && return local_tensor
  prev_vertex = only(prev_vertex_set)

  psi = state(problem)
  A = psi[prev_vertex]

  next_vertices = filter(v -> (it.hascommoninds(psi[v], A)), region)
  isempty(next_vertices) && return local_tensor
  next_vertex = only(next_vertices)
  C = psi[next_vertex]

  # Analyze indices of A
  # TODO: if "a" is missing, could supply a 1-dim index and put on both A and C?
  a = commonind(A, C)
  isnothing(a) && return local_tensor
  basis_inds = uniqueinds(A, C)

  ci = combinedind(combiner(basis_inds...))
  ax_space = expand_space(space(ci), expansion_factor)
  ax = Index(ax_space, "ax")

  linear_map(w) = (w - A * (dag(A) * w))
  Y = linear_map(random_itensor(basis_inds, dag(ax)))
  maxdim_goal = ceil(Int, expansion_factor * dim(a))
  maxdim_goal = min(maxdim_goal, maxdim)
  Ux_maxdim = maxdim_goal - dim(a)
  (norm(Y) <= 1E-15 || Ux_maxdim <= 0) && return local_tensor
  Ux, S, V = svd(Y, basis_inds; cutoff=1E-14, maxdim=Ux_maxdim, lefttags="ux,Link")

  Ux = linear_map(Ux)
  ux = commonind(Ux, S)
  Ax, sa = directsum(A => a, Ux => ux)
  expander = dag(Ax) * A
  psi[prev_vertex] = Ax
  psi[next_vertex] = expander * C
  local_tensor = expander*local_tensor

  return local_tensor
end
