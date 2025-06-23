using NDTensors: NDTensors
using NDTensors.BackendSelection: Backend, @Backend_str
import ConstructionBase: setproperties

default_expansion_factor() = 1.5
default_max_expand() = 4

function subspace_expand(
  problem,
  local_state,
  region_iterator;
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
  return subspace_expand(
    Backend(algorithm),
    problem,
    local_state,
    region_iterator;
    cutoff,
    mindim,
    maxdim,
    kws...,
  )
end

function subspace_expand(backend, problem, local_state, region_iterator; kws...)
  error(
    "Subspace expansion (subspace_expand!) not defined for requested combination of algorithm and problem types",
  )
end

function subspace_expand(
  backend::Backend{:nothing}, problem, local_state, region_iterator; kws...
)
  problem, local_state
end

function compute_expansion(
  current_dim,
  basis_size;
  expansion_factor=default_expansion_factor(),
  max_expand=default_max_expand(),
  maxdim=default_maxdim(),
)
  # Note: expand_maxdim will be *added* to current bond dimension
  # Obtain expand_maxdim from expansion_factor
  expand_maxdim = ceil(Int, expansion_factor * current_dim)
  # Enforce max_expand keyword
  expand_maxdim = min(max_expand, expand_maxdim)

  # Restrict expand_maxdim below theoretical upper limit
  expand_maxdim = min(basis_size-current_dim, expand_maxdim)
  # Enforce total maxdim setting (e.g. used in inserter step)
  expand_maxdim = min(maxdim-current_dim, expand_maxdim)
  # Ensure expand_maxdim is non-negative
  expand_maxdim = max(0, expand_maxdim)
  return expand_maxdim
end
