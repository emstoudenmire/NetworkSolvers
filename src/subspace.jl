using ITensors:
  commonind,
  dag,
  dim,
  directsum,
  dot,
  flux,
  hascommoninds,
  Index,
  inds,
  norm,
  onehot,
  pause,
  space,
  uniqueinds,
  random_itensor,
  qr

expand_space(χ::Integer, expansion_factor) = max(χ+1,floor(Int,expansion_factor*χ))

function expand_space(χs::Vector{<:Pair}, expansion_factor)
  return [q=>expand_space(d,expansion_factor) for (q,d) in χs]
end

total_space(spaces::Vector{<:Integer}) = prod(spaces)

function total_space(spaces::Vector)
  merge_space = Dict{QN,Int}()
  for it in Iterators.product(spaces)
    qtot = QN()
    dtot = 1
    for qd in it
      qtot += qd[1]
      dtot *= qd[2]
    end
    merge_space[qtot] += dtot
  end
  @show merge_space
  return [q=>d for (q,d) in merge_space]
end

function expand_space(a, basis_inds, expansion_factor)
  basis_space = total_space(space.(basis_inds))
  @show basis_space
  max_expand = dim_basis - dim(a)
end

function subspace_expand!(
  problem::EigsolveProblem, local_tensor, region; prev_region, expansion_factor=1.1, kws...
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

  #------------------------

  # Determine maximum value of num_expand
  #dim_basis = prod(dim.(basis_inds))
  #max_expand = dim_basis - dim(a)

  ax_space = expand_space(a,basis_inds,expansion_factor)
  pause()
  ax = Index(ax_space,"ax")
  #(num_expand <= 0) && return local_tensor

  linear_map(w) = (w - A * (dag(A) * w))

  Y = linear_map(random_itensor(basis_inds,ax))
  Q,R = qr(Y,basis_inds)
  q = commonind(Q,R)
  Ax, sa = directsum(A => a, Q => q)

  expander = dag(Ax) * A
  psi[prev_vertex] = Ax
  psi[next_vertex] = expander * C

  # TODO: avoid computing local tensor twice
  #       while also handling AbstractEdge region case
  local_tensor = prod(psi[v] for v in region)

  return local_tensor
end
