using NDTensors: NDTensors
using NDTensors.BackendSelection: Backend, @Backend_str

subspace_expand!(backend, problem, local_tensor, region; kws...) = local_tensor

function subspace_expand!(
  problem::EigsolveProblem,
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
