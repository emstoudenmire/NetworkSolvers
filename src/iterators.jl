import Base.Iterators: repeated

#
# sweep_iterator
#

function sweep_iterator(problem, sweep_kwargs_array)
  return [region_iterator(problem; sweep_kwargs...) for sweep_kwargs in sweep_kwargs_array]
end

sweep_iterator(problem, nsweeps::Integer) = sweep_iterator(problem,repeated((;),nsweeps))

#
# step_iterator
#

step_iterator(args...; kws...) = Iterators.flatten(sweep_iterator(args...; kws...))


#
# RegionIterator
#

@kwdef mutable struct RegionIterator{Problem,RegionPlan}
  problem::Problem
  region_plan::RegionPlan = region_plan(problem)
  which_region::Int = 1
  #extra_kwargs::NamedTuple = (;)
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

  R.problem = region_iterator_callback(problem(R); region=current_region, region_kwargs...)

  return R, next
end

function region_iterator(problem; sweep_kwargs...)
  return RegionIterator(; problem, region_plan=region_plan(problem; sweep_kwargs...))
end


function region_iterator_callback(problem; region, extracter_kwargs=(;), updater_kwargs=(;), inserter_kwargs=(;), kwargs...)
  problem, local_tensor = extracter(problem; region, kwargs...)
  problem, local_tensor = updater(problem, local_tensor; region, kwargs...)
  problem = inserter(problem, local_tensor, region; kwargs...)
  return problem
end

region_plan(problem; sweep_kwargs...) = basic_path_regions(state(problem); sweep_kwargs...)
