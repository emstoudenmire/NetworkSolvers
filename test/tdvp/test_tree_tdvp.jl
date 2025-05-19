import NetworkSolvers as ns
using Test: @test, @testset

using ITensors
import ITensorNetworks as itn
import Graphs as gr
import NamedGraphs as ng
import ITensorMPS as itm

function chain_plus_ancilla(; nchain)
  g = ng.NamedGraph()
  for j in 1:nchain
    gr.add_vertex!(g, j)
  end
  for j in 1:(nchain - 1)
    gr.add_edge!(g, j=>j+1)
  end
  # Add ancilla vertex near middle of chain
  gr.add_vertex!(g, 0)
  gr.add_edge!(g, 0=>nchain÷2)
  return g
end

@testset "Regression Test: Tree TDVP on chain plus ancilla" begin
  outputlevel = 1

  N = 10
  g = chain_plus_ancilla(; nchain=N)

  sites = itn.siteinds("S=1/2", g)

  # Make Heisenberg model Hamiltonian
  h = itm.OpSum()
  for j in 1:(N - 1)
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

  cutoff = 1E-10
  maxdim = 100
  nsweeps = 5

  nsites = 2
  E, gs_psi = ns.dmrg(
    H, psi0; inserter_kwargs=(; cutoff, maxdim), nsites, nsweeps, outputlevel
  )
  (outputlevel >= 1) && println("2-site DMRG energy = ", E)
  @show itn.maxlinkdim(gs_psi)

  inserter_kwargs=(; cutoff=1E-12, maxdim)
  nsites = 1
  tmax = 0.10
  time_range = 0.0:0.02:tmax
  psi1_t = ns.tdvp(H, gs_psi, time_range; inserter_kwargs, nsites, outputlevel=0)
  (outputlevel >= 1) && println("Done with $nsites-site TDVP")

  @test itn.norm(psi1_t) ≈ 1.0

  nsites = 2
  psi2_t = ns.tdvp(H, gs_psi, time_range; inserter_kwargs, nsites, outputlevel=0)
  (outputlevel >= 1) && println("Done with $nsites-site TDVP")
  @test itn.norm(psi2_t) ≈ 1.0

  @test abs(itn.inner(psi1_t, gs_psi)) > 0.99
  @test abs(itn.inner(psi1_t, psi2_t)) > 0.99

  # Test that accumulated phase angle is E*tmax
  z = itn.inner(psi1_t, gs_psi)
  @test atan(imag(z)/real(z)) ≈ E*tmax
end
