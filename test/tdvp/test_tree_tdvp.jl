import NetworkSolvers as ns
using Test: @test, @testset

using ITensors
import ITensorNetworks as itn
import Graphs as gr
import NamedGraphs as ng
import ITensorMPS as itm

function chain_plus_ancilla(; nchain)
  g = ng.NamedGraph()
  for j=1:nchain
    gr.add_vertex!(g, j)
  end
  for j=1:nchain-1
    gr.add_edge!(g, j=>j+1)
  end
  # Add ancilla vertex near middle of chain
  gr.add_vertex!(g, 0)
  gr.add_edge!(g, 0=>nchain÷2)
  return g
end

@testset "Regression Test: Tree TDVP on chain plus ancilla" begin
  outputlevel = 1

  N = 8
  g = chain_plus_ancilla(; nchain=N)

  sites = itn.siteinds("S=1/2", g)

  # Make Heisenberg model Hamiltonian
  h = itm.OpSum()
  for j in 1:N-1
    h += "Sz", j, "Sz", j+1
    h += 1/2, "S+", j, "S-", j+1
    h += 1/2, "S-", j, "S+", j+1
  end
  H = itn.ttn(h, sites)

  # Make initial product state
  state = Dict{Int,String}()
  for (j, v) in enumerate(gr.vertices(sites))
    state[v] = iseven(j) ? "Up" : "Dn"
  end
  psi0 = itn.ttn(state, sites)

  cutoff = 1E-5
  maxdim = 40
  nsweeps = 5

  nsites = 2
  E, gs_psi = ns.dmrg(
    H,
    psi0;
    truncation_kwargs=(; cutoff, maxdim),
    nsites,
    nsweeps,
    outputlevel,
  )
  (outputlevel >= 1) && println("2-site DMRG energy = ", E)


  nsites = 1
  time_range = 0.0:0.1:1.0
  psi_t = ns.tdvp(
    H,
    psi0,
    time_range;
    truncation_kwargs=(; cutoff, maxdim),
    nsites,
    outputlevel=0,
  )
  (outputlevel >= 1) && println("Done with $nsites-site TDVP")

end
