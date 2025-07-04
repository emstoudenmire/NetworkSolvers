module NetworkSolvers

using TensorOperations: TensorOperations # <-- to avoid a warning when
# loading ITensorNetworks

import ITensors as it
import ITensorNetworks as itn
import NamedGraphs as ng
using Graphs: Graphs

# Core solvers and linear algebra
include("solvers/eigsolve.jl")
include("solvers/exponentiate.jl")
include("solvers/runge_kutta.jl")
include("sketched_linear_algebra/range_finder.jl")
include("operator_map.jl")

# Algorithms and interfaces
include("eigsolve.jl")
include("applyexp.jl")
include("fitting.jl")

# Tensor network algorithm components
include("truncation_parameters.jl")
include("extracter.jl")
include("inserter.jl")
include("subspace/subspace.jl")
include("subspace/ortho_subspace.jl")
include("subspace/densitymatrix.jl")
include("permute_indices.jl")

# Iterators, tree traversal, region plans
include("iterators.jl")
include("adapters.jl")
include("region_plans/euler_tour.jl")
include("region_plans/euler_plans.jl")
include("region_plans/dfs_plans.jl")
include("region_plans/tdvp_region_plans.jl")
include("sweep_solve.jl")

end
