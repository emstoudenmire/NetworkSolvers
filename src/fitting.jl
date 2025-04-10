import ITensorNetworks as itn
import ITensors as it
import NamedGraphs.PartitionedGraphs as npg
using NamedGraphs: NamedEdge
using Printf

using ITensorNetworks:
  ITensorNetwork, induced_subgraph, ket_vertices, rename_vertices, inv_vertex_map

@kwdef mutable struct FittingProblem{State,OverlapNetwork}
  state::State
  overlapnetwork::OverlapNetwork
  overlap::Number = 0
  gauge_region
end

overlap(F::FittingProblem) = F.overlap
overlapnetwork(F::FittingProblem) = F.overlapnetwork
state(F::FittingProblem) = F.state
gauge_region(F::FittingProblem) = F.gauge_region

function extracter!(problem::FittingProblem, region; kws...)
  prev_region = gauge_region(problem)
  tn = state(problem)
  path = itn.edge_sequence_between_regions(tn, prev_region, region)
  tn = itn.gauge_walk(itn.Algorithm("orthogonalize"), tn, path)
  verts = unique(vcat(src.(path), dst.(path)))
  factors = [tn[v] for v in verts]
  o_tn = overlapnetwork(problem)
  o_tn = itn.update_factors(o_tn, Dict(zip(verts, factors)))
  pe_path = npg.partitionedges(itn.partitioned_tensornetwork(o_tn), path)
  o_tn = itn.update(
    itn.Algorithm("bp"), o_tn, pe_path; message_update_function_kwargs=(; normalize=false)
  )

  problem.state = tn
  problem.overlapnetwork = o_tn
  problem.gauge_region = region

  local_tensor = itn.environment(o_tn, region)
  sequence = itn.contraction_sequence(local_tensor; alg="optimal")
  local_tensor = dag(it.contract(local_tensor; sequence))

  return local_tensor
end

function prepare_subspace!(problem::FittingProblem, local_tensor, region; sweep, kws...)
  local_tensor = subspace_expand!(problem, local_tensor, region; sweep, kws...)
  return local_tensor
end

function inserter!(
  problem::FittingProblem,
  local_tensor,
  region;
  cutoff=default_cutoff(),
  maxdim=default_maxdim(),
  mindim=default_mindim(),
  normalize=true,
  sweep,
  kws...,
)
  cutoff = get_or_last(cutoff, sweep)
  mindim = get_or_last(mindim, sweep)
  maxdim = get_or_last(maxdim, sweep)

  psi = state(problem)
  v = last(region)
  if length(region) == 2
    e = ng.edgetype(psi)(first(region), last(region))
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
  normalize && (psi[v] /= norm(psi[v]))
  problem.state = psi
  #TODO: Why does this break?
  #problem.gauge_region = [v]

  on = overlapnetwork(problem)
  tn = state(problem)
  on = itn.update_factors(on, Dict(zip(region, [dag(tn[v]) for v in region])))
  problem.overlapnetwork = on

  return nothing
end

function updater!(F::FittingProblem, local_tensor, region; outputlevel, kws...)
  n = (local_tensor * dag(local_tensor))[]
  F.overlap = n / sqrt(n)
  if outputlevel >= 2
    @printf("  Region %s: squared overlap = %.12f\n", region, overlap(F))
  end
  return local_tensor
end

function fit_tensornetwork(
  overlap_network,
  args...;
  nsweeps=25,
  nsites=1,
  outputlevel=0,
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(; normalize=true),
  kws...,
)
  overlap_bpc = itn.BeliefPropagationCache(overlap_network, args...)
  ket_tn = first(induced_subgraph(overlap_network, ket_vertices(overlap_network)))
  init_prob = FittingProblem(;
    state=ket_tn, overlapnetwork=overlap_bpc, gauge_region=collect(vertices(ket_tn))
  )

  common_sweep_kwargs = (; nsites, outputlevel, updater_kwargs, inserter_kwargs)
  kwargs_array = [(; common_sweep_kwargs..., sweep=s) for s in 1:nsweeps]
  sweep_iter = sweep_iterator(init_prob, kwargs_array)
  converged_prob = alternating_update(sweep_iter; outputlevel, kws...)
  return rename_vertices(inv_vertex_map(overlap_network), state(converged_prob))
end

function fit_tensornetwork(tn, init_state, args...; kwargs...)
  return fit_tensornetwork(itn.inner_network(tn, init_state), args; kwargs...)
end

function itn.truncate(tn; maxdim::Int64, cutoff=0.0, kwargs...)
  init_state = itn.ITensorNetwork(
    v -> inds -> it.delta(inds), itn.siteinds(tn); link_space=maxdim
  )
  overlap_network = itn.inner_network(tn, init_state)
  return fit_tensornetwork(overlap_network; inserter_kwargs=(; cutoff, maxdim), kwargs...)
end

function itn.apply(
  A::ITensorNetwork, x::ITensorNetwork; maxdim::Int64, cutoff=0.0, kwargs...
)
  init_state = itn.ITensorNetwork(
    v -> inds -> it.delta(inds), itn.siteinds(x); link_space=maxdim
  )
  overlap_network = itn.inner_network(x, A, init_state)
  return fit_tensornetwork(overlap_network; inserter_kwargs=(; cutoff, maxdim), kwargs...)
end
