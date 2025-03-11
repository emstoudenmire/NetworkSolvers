import ITensors as it
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
import NetworkSolvers as ns

function main(; N=10)
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
      @printf("  Region %s: energy = %.12f\n",region,data.current_energy)
    end
    println("Done with sweep $sweep")
  end

  return nothing
end
