import ITensors as it
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
#import NetworkSolvers as ns
using Random
using Graphs: edges, dst, src, vertices
using NamedGraphs.NamedGraphGenerators: named_comb_tree

include("timed_eigsolve.jl")

function main(; N=100, nsites=2, site_type="S=1", conserve_qns=false)
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

  nsweeps = 4
  cutoff = 1E-9
  maxdim = [10, 40, 80, 160]
  outputlevel = 1
  truncation_kwargs = (; cutoff, maxdim)
  subspace_kwargs = (;) #(; algorithm="densitymatrix", max_expand=4)

  @time begin
    energy, gs_psi = timed_dmrg(
      H, psi; nsweeps, nsites, truncation_kwargs, subspace_kwargs, outputlevel
    )
  end
  println("Final energy = ", energy)

  if site_type == "S=1" && N == 10
    println("Exact energy = -12.8945601")
  end

  return nothing
end
