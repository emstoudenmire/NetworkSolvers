import ITensorNetworks as itn
import ITensors as it
import Graphs: vertices
import NamedGraphs.PartitionedGraphs as npg
using NamedGraphs: NamedEdge
using Printf

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

function set!(
  F::FittingProblem;
  state=state(F),
  overlapnetwork=overlapnetwork(F),
  gauge_region=gauge_region(F),
  overlap=overlap(F),
)
  F = FittingProblem(; state, overlapnetwork, gauge_region, overlap)
end

function extracter!(problem::FittingProblem, region_iterator; kws...)
  region = current_region(region_iterator)
  prev_region = gauge_region(problem)
  tn = state(problem)
  o_tn = itn.update_factors(
    overlapnetwork(problem), Dict(zip(region, [tn[v] for v in prev_region]))
  )
  path = itn.edge_sequence_between_regions(tn, prev_region, region)
  tn = itn.gauge_walk(itn.Algorithm("orthogonalize"), tn, path)
  verts = unique(vcat(src.(path), dst.(path)))
  factors = [tn[v] for v in verts]
  o_tn = itn.update_factors(o_tn, Dict(zip(verts, factors)))
  pe_path = npg.partitionedges(itn.partitioned_tensornetwork(o_tn), path)
  o_tn = itn.update(
    itn.Algorithm("bp"), o_tn, pe_path; message_update_function_kwargs=(; normalize=false)
  )
  set!(problem; state=tn, overlapnetwork=o_tn, gauge_region=region)

  local_tensor = itn.environment(o_tn, region)
  sequence = itn.contraction_sequence(local_tensor; alg="optimal")
  local_tensor = dag(it.contract(local_tensor; sequence))

  return local_tensor
end

function prepare_subspace!(problem::FittingProblem, local_tensor, region; sweep, kws...)
  local_tensor = subspace_expand!(problem, local_tensor, region; sweep, kws...)
  return local_tensor
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
  truncation_kwargs=(;),
  normalize=true,
  kws...,
)
  overlap_bpc = itn.BeliefPropagationCache(overlap_network, args...)
  ket_tn = first(itn.induced_subgraph(overlap_network, itn.ket_vertices(overlap_network)))
  init_prob = FittingProblem(;
    state=ket_tn, overlapnetwork=overlap_bpc, gauge_region=collect(vertices(ket_tn))
  )

  truncation_kwargs = (; truncation_kwargs..., normalize, set_orthogonal_region=false)
  common_sweep_kwargs = (; nsites, outputlevel, updater_kwargs, truncation_kwargs)
  kwargs_array = [(; common_sweep_kwargs..., sweep=s) for s in 1:nsweeps]
  sweep_iter = sweep_iterator(init_prob, kwargs_array)
  converged_prob = alternating_update(sweep_iter; outputlevel, kws...)
  return itn.rename_vertices(itn.inv_vertex_map(overlap_network), state(converged_prob))
end

function fit_tensornetwork(tn, init_state, args...; kwargs...)
  return fit_tensornetwork(itn.inner_network(tn, init_state), args; kwargs...)
end

function itn.truncate(tn; maxdim::Int64, cutoff=0.0, kwargs...)
  init_state = itn.ITensorNetwork(
    v -> inds -> it.delta(inds), itn.siteinds(tn); link_space=maxdim
  )
  overlap_network = itn.inner_network(tn, init_state)
  return fit_tensornetwork(overlap_network; truncation_kwargs=(; cutoff, maxdim), kwargs...)
end

function itn.apply(
  A::itn.ITensorNetwork, x::itn.ITensorNetwork; maxdim::Int64, cutoff=0.0, kwargs...
)
  init_state = itn.ITensorNetwork(
    v -> inds -> it.delta(inds), itn.siteinds(x); link_space=maxdim
  )
  overlap_network = itn.inner_network(x, A, init_state)
  return fit_tensornetwork(overlap_network; truncation_kwargs=(; cutoff, maxdim), kwargs...)
end
