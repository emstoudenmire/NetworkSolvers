import ITensorNetworks as itn
import NetworkSolvers as ns
using NamedGraphs.NamedGraphGenerators: named_comb_tree
import NetworkSolvers as ns
using Test: @test, @testset
using Printf
using StableRNGs: StableRNG
using ITensorNetworks.ModelHamiltonians: heisenberg
import ITensorMPS as itm

# This part of the code is just running DMRG to obtain a 
# reasonable state to apply an operator to
function compute_state(N, sites)
  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, sites)
  state = Dict{Int,String}()
  for v in itn.vertices(H)
    state[v] = iseven(v) ? "Up" : "Dn"
  end
  psi0 = itn.ttn(state, sites)
  trunc = (; maxdim=50, cutoff=1E-5)
  energy, psi = ns.dmrg(
    H,
    psi0;
    nsweeps=5,
    nsites=2,
    extracter_kwargs=(; trunc),
    inserter_kwargs=(; trunc),
    outputlevel=1,
  )
  return psi
end

@testset "Fitting Regression Test" begin
  N = 8
  conserve_qns = false
  sites = itn.siteinds("S=1/2", N; conserve_qns)

  terms = itm.OpSum()
  terms += "S+", 3, "S-", 5
  terms += "S-", 3, "S+", 5
  O = itn.ITensorNetwork(itn.ttn(terms, sites))

  # Issue #1, apply/fitting doesn't work for TreeTensorNetwork type
  psi = compute_state(N, sites)
  psi = itn.ITensorNetwork(psi)

  # Issue #2, apply/fitting doesn't work with QNs
  Opsi = itn.apply(O, psi; maxdim=60, nsites=2, normalize=false)
  f = itn.inner(Opsi, psi; alg="exact") / itn.inner(psi, O, psi; alg="exact")
  @printf("One-site apply. Fidelity = %s\n", f)
end
