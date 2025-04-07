import ITensorNetworks as itn
import ITensors as it
import NamedGraphs.PartitionedGraphs as npg
using NamedGraphs: NamedEdge
using Printf

using ITensorNetworks: ITensorNetwork

@kwdef mutable struct FittingProblem{State,OverlapNetwork}
  state::State
  overlapnetwork::OverlapNetwork
  squared_scalar::Number = 0
end

squared_scalar(F::FittingProblem) = F.squared_scalar
overlapnetwork(F::FittingProblem) = F.overlapnetwork
state(F::FittingProblem) = F.state

function set(
  F::FittingProblem;
  state=state(F),
  overlapnetwork=overlapnetwork(F),
  squared_scalar=squared_scalar(F),
)
  return FittingProblem(; state, overlapnetwork, squared_scalar)
end

function extracter!(problem::FittingProblem, prev_region, region; kws...)
  tn = state(problem)
  path = itn.edge_sequence_between_regions(tn, prev_region, region)
  tn = itn.gauge_walk(itn.Algorithm("orthogonalize"), tn, path)
  verts = unique(vcat(src.(path), dst.(path)))
  factors = [tn[v] for v in verts]
  o_tn = overlapnetwork(problem)
  o_tn = itn.update_factors(o_tn, Dict(zip([(v, "ket") for v in verts], factors)))
  pe_path = npg.partitionedges(
    itn.partitioned_tensornetwork(o_tn),
    [NamedEdge((src(e), "ket") => (dst(e), "ket")) for e in path],
  )
  o_tn = itn.update(
    itn.Algorithm("bp"), o_tn, pe_path; message_update_function_kwargs=(; normalize=false)
  )

  problem.state = tn
  problem.overlapnetwork = o_tn

  local_tensor = itn.environment(o_tn, [(v, "ket") for v in region])
  sequence = itn.contraction_sequence(local_tensor; alg="optimal")
  local_tensor = dag(it.contract(local_tensor; sequence))

  return local_tensor
end

function fitting_inserter!(problem::FittingProblem, local_tensor, region; kws...)
  inserter!(problem, local_tensor, region; kws...)
  on = overlapnetwork(problem)
  tn = state(problem)
  on = itn.update_factors(
    on, Dict(zip([(v, "ket") for v in region], [dag(tn[v]) for v in region]))
  )
  problem.overlapnetwork = on
  return nothing
end

function region_iterator_action!(
  problem::FittingProblem;
  region,
  prev_region=nothing,
  extracter_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  sweep,
  kwargs...,
)
  prev_region = isnothing(prev_region) ? collect(vertices(state(problem))) : prev_region
  local_tensor = extracter!(problem, prev_region, region; extracter_kwargs..., kwargs...)
  updater!(problem, local_tensor, region; updater_kwargs..., kwargs...)
  fitting_inserter!(
    problem,
    local_tensor,
    region;
    set_orthogonal_region=false,
    sweep,
    inserter_kwargs...,
    kwargs...,
  )
  return nothing
end

function updater!(F::FittingProblem, local_tensor, region; outputlevel, kws...)
  n = (local_tensor * dag(local_tensor))[]
  F.squared_scalar = n / sqrt(n)
  if outputlevel >= 2
    @printf("  Region %s: squared overlap = %.12f\n", region, squared_scalar(F))
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
  init_prob = FittingProblem(;
    state=itn.ket_network(overlap_network), overlapnetwork=overlap_bpc
  )
  common_sweep_kwargs = (; nsites, outputlevel, updater_kwargs, inserter_kwargs)
  kwargs_array = [(; common_sweep_kwargs..., sweep=s) for s in 1:nsweeps]
  sweep_iter = sweep_iterator(init_prob, kwargs_array)
  converged_prob = alternating_update(sweep_iter; outputlevel, kws...)
  return state(converged_prob)
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
