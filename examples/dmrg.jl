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

  nsweeps = 5
  cutoff = 1E-12
  maxdim = [10, 40, 80, 160]
  outputlevel = 1
  trunc = (; cutoff, maxdim)
  extracter_kwargs = (; trunc, subspace_algorithm="densitymatrix", expansion_factor=1.1)
  inserter_kwargs = (; trunc)

  @time begin
    energy, gs_psi = ns.dmrg(
      H, psi; nsweeps, nsites, extracter_kwargs, inserter_kwargs, outputlevel
    )
  end
  println("Final energy = ", energy)

  if site_type == "S=1" && N == 10
    println("Exact energy = -12.8945601")
  end

  return nothing
end

function tree_dmrg(; conserve_qns=false)
  tooth_lengths = [5, 5, 5]
  c = named_comb_tree(tooth_lengths)
  s = itn.siteinds("S=1", c; conserve_qns)

  os = itm.OpSum()
  for e in edges(c)
    os += "Sz", src(e), "Sz", dst(e)
    os += 1 / 2, "S+", src(e), "S-", dst(e)
    os += 1 / 2, "S-", src(e), "S+", dst(e)
  end
  H = itn.ttn(os, s)

  psi = itn.random_mps(s; link_space=4)

  nsweeps = 14
  nsites = 2
  cutoff = 1E-9
  maxdim = 10
  outputlevel = 2
  trunc = (; cutoff, maxdim)
  extracter_kwargs = (; trunc, subspace_algorithm="densitymatrix", expansion_factor=1.2)
  inserter_kwargs = (; trunc)
  energy, gs_psi = ns.dmrg(
    H, psi; nsweeps, nsites, extracter_kwargs, inserter_kwargs, outputlevel
  )
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
