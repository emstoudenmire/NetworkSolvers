using ITensors
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
using Random
import NetworkSolvers as ns

function main(; N=10, total_time=1.0, time_step=0.1)
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi = itn.random_mps(s; link_space=4)

  cutoff = 1E-6
  maxdim = [10, 20, 40, 100, 200]
  outputlevel = 2
  inserter_kwargs = (; maxdim, cutoff)
  res = ns.applyexp(H, psi, 0.0:time_step:total_time; inserter_kwargs, outputlevel)

  return nothing
end

function test_tdvp()
  Random.seed!(1)
  total_time = 1.0
  time_step = 0.10

  N = 6
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi = itn.random_mps(s; link_space=32)

  outputlevel = 1

  nsweeps = convert(Int, ceil(abs(total_time / time_step)))
  @assert (nsweeps * time_step ≈ total_time)
  tdvp_iter = tdvp_region_iterator(H, psi; time_step)

  szs = zeros(nsweeps, N)
  for sweep in 1:nsweeps
    for i in tdvp_iter
    end
    szs[sweep, :] = itn.expect(data(tdvp_iter).state, "Sz")
  end
  println("\nResult from TDVP:")
  display(szs)

  #
  # Use ED to check
  #
  H = prod(H)
  psi = prod(psi)
  expH = exp(H * time_step)
  szs_ed = zeros(nsweeps, N)
  for sweep in 1:nsweeps
    psi = noprime(psi * expH)
    psi /= norm(psi)
    for j in 1:N
      szs_ed[sweep, j] = scalar(dag(prime(psi, s[j])) * itm.op("Sz", s[j]) * psi)
    end
  end
  println("\nResult from ED:")
  display(szs_ed)

  return nothing
end
