using NDTensors: NDTensors
using NDTensors.BackendSelection: Backend, @Backend_str

default_expansion_factor() = 1.5
default_max_expand() = 4

function prepare_subspace!(
  problem, local_tensor, region; prev_region=nothing, sweep, kws...
)
  local_tensor = subspace_expand!(problem, local_tensor, region; prev_region, sweep, kws...)
  problem.operator = itn.position(operator(problem), state(problem), region)
  return local_tensor
end

subspace_expand!(backend, problem, local_tensor, region; kws...) = local_tensor

function subspace_expand!(
  problem,
  local_tensor,
  region;
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  algorithm=nothing,
  sweep,
  kws...,
)
  cutoff = get_or_last(cutoff, sweep)
  mindim = get_or_last(mindim, sweep)
  maxdim = get_or_last(maxdim, sweep)
  return subspace_expand!(
    Backend(algorithm), problem, local_tensor, region; cutoff, mindim, maxdim, kws...
  )
end

function compute_expansion(
  current_dim,
  basis_size;
  expansion_factor=default_expansion_factor(),
  max_expand=default_max_expand(),
  maxdim=default_maxdim(),
)
  expand_maxdim = ceil(Int, expansion_factor * current_dim)
  expand_maxdim = min(basis_size-current_dim, expand_maxdim)
  expand_maxdim = min(maxdim-current_dim, expand_maxdim)
  expand_maxdim = max(0, expand_maxdim)
  return expand_maxdim
end
