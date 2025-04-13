import ITensors as it
import ITensorNetworks as itn
using Printf
import NetworkSolvers as ns
using Random
using Graphs: edges, dst, src
using NamedGraphs.NamedGraphGenerators: named_grid, named_comb_tree
using TensorOperations: TensorOperations
using StableRNGs: StableRNG
using ITensorNetworks: ITensorNetwork, apply, inner, siteinds, random_tensornetwork, ttn
using ITensorNetworks.ModelHamiltonians: heisenberg

it.disable_warn_order()

function main()
  for elt in (Float32, Float64, Complex{Float32}, Complex{Float64})
    println();
    @show elt
    g = named_comb_tree((3, 2))
    s = siteinds("S=1/2", g)

    rng = StableRNG(1234)
    #One-site truncation
    a = random_tensornetwork(rng, elt, s; link_space=3)
    b = truncate(a; maxdim=3)
    f = inner(a, b; alg="exact") / sqrt(inner(a, a; alg="exact") * inner(b, b; alg="exact"))
    @printf("One-site truncation. Fidelity = %s\n", f)
    @assert abs(abs(f) - 1.0) <= 10*eps(real(elt))

    #Two-site truncation
    a = random_tensornetwork(rng, elt, s; link_space=3)
    b = truncate(a; maxdim=3, cutoff=1e-16, nsites=2)
    f = inner(a, b; alg="exact") / sqrt(inner(a, a; alg="exact") * inner(b, b; alg="exact"))
    @printf("Two-site truncation. Fidelity = %s\n", f)
    @assert abs(abs(f) - 1.0) <= 10*eps(real(elt))

    # #One-site apply (no normalization)
    a = random_tensornetwork(rng, elt, s; link_space=2)
    H = ITensorNetwork(ttn(heisenberg(g), s))
    Ha = apply(H, a; maxdim=4, nsites=1, normalize=false)
    f = inner(Ha, a; alg="exact") / inner(a, H, a; alg="exact")
    @printf("One-site apply. Fidelity = %s\n", f)
    @assert abs(f - 1.0) <= 10*eps(real(elt))

    # #Two-site apply (no normalization)
    a = random_tensornetwork(rng, elt, s; link_space=2)
    H = ITensorNetwork(ttn(heisenberg(g), s))
    Ha = apply(H, a; maxdim=4, cutoff=1e-16, nsites=2, normalize=false)
    f = inner(Ha, a; alg="exact") / inner(a, H, a; alg="exact")
    @printf("Two-site apply. Fidelity = %s\n", f)
    @assert abs(f - 1.0) <= 10*eps(real(elt))
  end
end

main()
