import NetworkSolvers as ns
using Test: @test, @testset

using ITensors
import ITensorNetworks as itn
import Graphs as gr
import NamedGraphs as ng
import ITensorMPS as itm

include("../utilities/simple_ed_methods.jl")

"""
  build_tree

  Make a tree with central vertex (0,0) and
  nbranch branches of nbranch_sites each.
"""
function build_tree(; nbranch=3, nbranch_sites=3)
  g = ng.NamedGraph()
  gr.add_vertex!(g, (0, 0))
  for branch in 1:nbranch, site in 1:nbranch_sites
    gr.add_vertex!(g, (branch, site))
  end
  for branch in 1:nbranch
    gr.add_edge!(g, (0, 0)=>(branch, 1))
    for site in 2:nbranch_sites
      gr.add_edge!(g, (branch, site-1)=>(branch, site))
    end
  end
  return g
end

@testset "Tree DMRG" begin
  outputlevel = 1

  g = build_tree(; nbranch=3, nbranch_sites=3)

  sites = itn.siteinds("S=1/2", g)

  # Make Heisenberg model Hamiltonian
  h = itm.OpSum()
  for edge in gr.edges(sites)
    i, j = gr.src(edge), gr.dst(edge)
    h += "Sz", i, "Sz", j
    h += 1/2, "S+", i, "S-", j
    h += 1/2, "S-", i, "S+", j
  end
  H = itn.ttn(h, sites)

  # Make initial product state
  state = Dict{Tuple{Int,Int},String}()
  for (j, v) in enumerate(gr.vertices(sites))
    state[v] = iseven(j) ? "Up" : "Dn"
  end
  psi0 = itn.ttn(state, sites)

  (outputlevel >= 1) && println("Computing exact ground state")
  Ex, psix = ed_ground_state(H,psi0)
  (outputlevel >= 1) && println("Ex = ",Ex)

  cutoff = 1E-5
  maxdim = 40
  nsweeps = 5

  #
  # Test 2-site DMRG without subspace expansion
  #
  nsites = 2
  subspace_kwargs = (;)#(; algorithm="densitymatrix")
  E, psi = ns.dmrg(
    H,
    psi0;
    inserter_kwargs=(; cutoff, maxdim),
    nsites,
    nsweeps,
    subspace_kwargs,
    outputlevel,
  )
  (outputlevel >= 1) && println("2-site DMRG energy = ",E)
  @test abs(E-Ex) < 1E-5


  #
  # Test 1-site DMRG with subspace expansion
  #
  nsites = 1
  nsweeps = 10
  subspace_kwargs = (; algorithm="densitymatrix", max_expand=8)
  cutoff = 1E-10
  maxdim = 200
  E, psi = ns.dmrg(
    H,
    psi0;
    inserter_kwargs=(; cutoff, maxdim),
    nsites,
    nsweeps,
    subspace_kwargs,
    outputlevel,
  )
  (outputlevel >= 1) && println("1-site+subspace DMRG energy = ",E)
  @test abs(E-Ex) < 1E-5

end
