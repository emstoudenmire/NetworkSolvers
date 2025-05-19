import ITensorNetworks as itn
import NetworkSolvers as ns
using NamedGraphs.NamedGraphGenerators: named_comb_tree
import NetworkSolvers as ns
using Test: @test, @testset
using Printf
using StableRNGs: StableRNG
using ITensorNetworks.ModelHamiltonians: heisenberg

@testset "Fitting Tests" begin
  outputlevel = 1
  for elt in (Float32, Float64, Complex{Float32}, Complex{Float64})
    (outputlevel >= 1) && println("\nFitting tests with elt = ", elt)
    g = named_comb_tree((3, 2))
    s = itn.siteinds("S=1/2", g)

    rng = StableRNG(1234)
    #One-site truncation
    a = itn.random_tensornetwork(rng, elt, s; link_space=3)
    b = truncate(a; maxdim=3)
    f =
      itn.inner(a, b; alg="exact") /
      sqrt(itn.inner(a, a; alg="exact") * itn.inner(b, b; alg="exact"))
    (outputlevel >= 1) && @printf("One-site truncation. Fidelity = %s\n", f)
    @test abs(abs(f) - 1.0) <= 10*eps(real(elt))

    #Two-site truncation
    a = itn.random_tensornetwork(rng, elt, s; link_space=3)
    b = truncate(a; maxdim=3, cutoff=1e-16, nsites=2)
    f =
      itn.inner(a, b; alg="exact") /
      sqrt(itn.inner(a, a; alg="exact") * itn.inner(b, b; alg="exact"))
    (outputlevel >= 1) && @printf("Two-site truncation. Fidelity = %s\n", f)
    @test abs(abs(f) - 1.0) <= 10*eps(real(elt))

    # #One-site apply (no normalization)
    a = itn.random_tensornetwork(rng, elt, s; link_space=2)
    H = itn.ITensorNetwork(itn.ttn(heisenberg(g), s))
    Ha = itn.apply(H, a; maxdim=4, nsites=1, normalize=false)
    f = itn.inner(Ha, a; alg="exact") / itn.inner(a, H, a; alg="exact")
    (outputlevel >= 1) && @printf("One-site apply. Fidelity = %s\n", f)
    @test abs(f - 1.0) <= 10*eps(real(elt))

    # #Two-site apply (no normalization)
    a = itn.random_tensornetwork(rng, elt, s; link_space=2)
    H = itn.ITensorNetwork(itn.ttn(heisenberg(g), s))
    Ha = itn.apply(H, a; maxdim=4, cutoff=1e-16, nsites=2, normalize=false)
    f = itn.inner(Ha, a; alg="exact") / itn.inner(a, H, a; alg="exact")
    (outputlevel >= 1) && @printf("Two-site apply. Fidelity = %s\n", f)
    @test abs(f - 1.0) <= 10*eps(real(elt))
  end
end
