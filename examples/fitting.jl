using ITensors
import ITensorNetworks as itn
using Printf
import NetworkSolvers as ns
using NamedGraphs.NamedGraphGenerators: named_comb_tree
using TensorOperations: TensorOperations
using StableRNGs: StableRNG
using ITensorNetworks.ModelHamiltonians: heisenberg

ITensors.disable_warn_order()

function fitting()
  for elt in (Float32, Float64, Complex{Float32}, Complex{Float64})
    println()
    @show elt
    g = named_comb_tree((3, 2))
    s = itn.siteinds("S=1/2", g)

    rng = StableRNG(1234)
    #One-site truncation
    a = itn.random_tensornetwork(rng, elt, s; link_space=3)
    b = ns.truncate(a; maxdim=3)
    f =
      itn.inner(a, b; alg="exact") /
      sqrt(itn.inner(a, a; alg="exact") * itn.inner(b, b; alg="exact"))
    @printf("One-site truncation. Fidelity = %s\n", f)
    @assert abs(abs(f) - 1.0) <= 10*eps(real(elt))

    #Two-site truncation
    a = itn.random_tensornetwork(rng, elt, s; link_space=3)
    b = ns.truncate(a; maxdim=3, cutoff=1e-16, nsites=2)
    f =
      itn.inner(a, b; alg="exact") /
      sqrt(itn.inner(a, a; alg="exact") * itn.inner(b, b; alg="exact"))
    @printf("Two-site truncation. Fidelity = %s\n", f)
    @assert abs(abs(f) - 1.0) <= 10*eps(real(elt))

    # #One-site apply (no normalization)
    a = itn.random_tensornetwork(rng, elt, s; link_space=2)
    H = itn.ITensorNetwork(itn.ttn(heisenberg(g), s))
    Ha = ns.apply(H, a; maxdim=4, nsites=1, normalize=false)
    f = itn.inner(Ha, a; alg="exact") / itn.inner(a, H, a; alg="exact")
    @printf("One-site apply. Fidelity = %s\n", f)
    @assert abs(f - 1.0) <= 10*eps(real(elt))

    # #Two-site apply (no normalization)
    a = itn.random_tensornetwork(rng, elt, s; link_space=2)
    H = itn.ITensorNetwork(itn.ttn(heisenberg(g), s))
    Ha = ns.apply(H, a; maxdim=4, cutoff=1e-16, nsites=2, normalize=false)
    f = itn.inner(Ha, a; alg="exact") / itn.inner(a, H, a; alg="exact")
    @printf("Two-site apply. Fidelity = %s\n", f)
    @assert abs(f - 1.0) <= 10*eps(real(elt))
  end
end
