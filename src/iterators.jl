
#
# sweep_iterator
#

function sweep_iterator(problem, sweep_kws)
  return [region_iterator(problem; sweep=s, kws...) for (s, kws) in enumerate(sweep_kws)]
end

function sweep_iterator(problem, nsweeps::Integer; sweep_kws...)
  return sweep_iterator(problem, Iterators.repeated(sweep_kws, nsweeps))
end

#
# step_iterator
#

step_iterator(args...; kws...) = Iterators.flatten(sweep_iterator(args...; kws...))

#
# RegionIterator
#

@kwdef mutable struct RegionIterator{Problem,RegionPlan}
  problem::Problem
  region_plan::RegionPlan
  which_region::Int = 1
end

problem(R::RegionIterator) = R.problem
current_region_plan(R::RegionIterator) = R.region_plan[R.which_region]
current_region(R::RegionIterator) = R.region_plan[R.which_region][1]
region_kwargs(R::RegionIterator) = R.region_plan[R.which_region][2]

function Base.iterate(R::RegionIterator, which=1)
  R.which_region = which
  region_plan_state = iterate(R.region_plan, which)
  isnothing(region_plan_state) && return nothing
  (current_region, region_kwargs), next = region_plan_state

  region_iterator_action!(problem(R); region=current_region, region_kwargs...)
  return R, next
end

#
# Functions associated with RegionIterator
#

function region_iterator(problem; nsites=1, sweep_kwargs...)
  return RegionIterator(;
    problem, region_plan=region_plan(problem; nsites, sweep_kwargs...)
  )
end

function region_iterator_action!(
  problem;
  region,
  extracter_kwargs=(;),
  subspace_kwargs=(;),
  updater_kwargs=(;),
  truncation_kwargs=(;),
  sweep,
  kwargs...,
)
  local_tensor = extracter!(problem, region; extracter_kwargs..., kwargs...)
  local_tensor = prepare_subspace!(
    problem, local_tensor, region; subspace_kwargs..., sweep, kwargs...
  )
  local_tensor = updater!(problem, local_tensor, region; updater_kwargs..., kwargs...)
  inserter!(problem, local_tensor, region; sweep, truncation_kwargs..., kwargs...)
  return nothing
end

function region_plan(problem; nsites, sweep_kwargs...)
  return basic_region_plan(state(problem); nsites, sweep_kwargs...)
end
