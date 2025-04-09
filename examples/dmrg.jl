import ITensors as it
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
import NetworkSolvers as ns
using Random
using Graphs: edges, dst, src, vertices
using NamedGraphs.NamedGraphGenerators: named_comb_tree

function dmrg(; N=10, nsites=2, site_type="S=1", conserve_qns=false)
  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  s = itn.siteinds(site_type, N; conserve_qns)
  H = itn.mpo(os, s)

  state = Dict{Int,String}()
  for v in vertices(H)
    state[v] = iseven(v) ? "Up" : "Dn"
  end
  psi = itn.ttn(state, s)

  nsweeps = 8
  cutoff = 1E-9
  maxdim = [10, 40, 80, 160]
  outputlevel = 1
  inserter_kwargs = (; cutoff, maxdim)
  subspace_kwargs = (; algorithm="densitymatrix", maxdim=4)

  @time begin
    energy, gs_psi = ns.dmrg(
      H, psi; nsweeps, nsites, inserter_kwargs, subspace_kwargs, outputlevel
    )
  end
  println("Final energy = ", energy)

  if site_type == "S=1" && N == 10
    println("Exact energy = -12.8945601")
  end

  return nothing
end

function tree_dmrg()
  tooth_lengths = [5, 5, 5]
  c = named_comb_tree(tooth_lengths)
  s = itn.siteinds("S=1", c)

  os = itm.OpSum()
  for e in edges(c)
    os += "Sz", src(e), "Sz", dst(e)
    os += 1 / 2, "S+", src(e), "S-", dst(e)
    os += 1 / 2, "S-", src(e), "S+", dst(e)
  end
  H = itn.ttn(os, s)

  psi = itn.random_mps(s; link_space=4)

  nsweeps = 14
  cutoff = 1E-9
  maxdim = 50
  outputlevel = 2
  nsites = 1
  inserter_kwargs = (; cutoff, maxdim)
  # TODO: 1-site and 2-site DMRG giving different energies
  #       is 1 site skipping a site? Is 2 site skipping any bonds?
  energy, gs_psi = ns.dmrg(H, psi; nsweeps, nsites, inserter_kwargs, outputlevel)
  return println("Final energy = ", energy)
end

function sweep_loop_version()
  N = 10
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi = itn.random_mps(s; link_space=4)

  nsweeps = 2
  cutoff = 1E-6
  maxdim = [10, 20, 40, 100, 200]
  outputlevel = 1
  region_iterator = dmrg_region_iterator(H, psi; nsweeps, maxdim, cutoff, outputlevel)

  for sweep in 1:nsweeps
    println("\nSweep $sweep:")
    for (region, data) in region_iterator
      @printf("  Region %s: energy = %.12f\n", region, data.current_energy)
    end
    println("Done with sweep $sweep")
  end

  return nothing
end
