
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
current_region(R::RegionIterator) = current_region_plan(R)[1]
region_kwargs(R::RegionIterator) = current_region_plan(R)[2]
function previous_region(R::RegionIterator)
  R.which_region==1 ? nothing : R.region_plan[R.which_region - 1][1]
end
function next_region(R::RegionIterator)
  R.which_region==length(R.region_plan) ? nothing : R.region_plan[R.which_region + 1][1]
end
is_last_region(R::RegionIterator) = isnothing(next_region(R))

function Base.iterate(R::RegionIterator, which=1)
  R.which_region = which
  region_plan_state = iterate(R.region_plan, which)
  isnothing(region_plan_state) && return nothing
  (current_region, region_kwargs), next = region_plan_state
  R.problem = region_iterator_action(problem(R), R; region_kwargs...)
  return R, next
end

#
# Functions associated with RegionIterator
#

function region_iterator(problem; sweep_kwargs...)
  return RegionIterator(; problem, region_plan=region_plan(problem; sweep_kwargs...))
end

function region_iterator_action(
  problem,
  region_iterator;
  extracter_kwargs=(;),
  subspace_kwargs=(;),
  updater_kwargs=(;),
  inserter_kwargs=(;),
  sweep,
  kwargs...,
)
  problem, local_state = extracter(
    problem, region_iterator; extracter_kwargs..., subspace_kwargs..., sweep, kwargs...
  )
  problem, local_state = updater(
    problem, local_state, region_iterator; updater_kwargs..., kwargs...
  )
  problem = inserter(
    problem, local_state, region_iterator; sweep, inserter_kwargs..., kwargs...
  )
  return problem
end

function region_plan(problem; nsites, sweep_kwargs...)
  return euler_sweep(state(problem); nsites, sweep_kwargs...)
end
