using ITensors
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
using Random
import NetworkSolvers as ns

function rk4_solver(H, t, ψ0; kws...)
  k1 = H(ψ0)
  k2 = k1 + (t / 2) * H(k1)
  k3 = k1 + (t / 2) * H(k2)
  k4 = k1 + t * H(k3)
  return ψ0 + (t / 6) * (k1 + 2 * k2 + 2 * k3 + k4), (;)
end

function main(; N=10, total_time=-1.0, time_step=-0.1)
  Random.seed!(1)
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi0 = itn.random_mps(s; link_space=16)

  cutoff = 1E-10
  maxdim = 16
  outputlevel = 0
  inserter_kwargs = (; maxdim, cutoff, normalize=true)
  time_range = 0.0:time_step:total_time
  res = ns.applyexp(H, psi0, time_range; inserter_kwargs, outputlevel)

  # Using RK solver
  updater_kwargs = (; solver=rk4_solver)
  res_rk4 = ns.applyexp(
    H, psi0, 0.0:time_step:total_time; inserter_kwargs, updater_kwargs, outputlevel
  )

  @show inner(res, res_rk4)

  return nothing
end

function test_tdvp(; N=6, total_time=-0.1, time_step=-0.01)
  Random.seed!(1)
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi = itn.random_mps(s; link_space=32)

  time_range = 0.0:time_step:total_time
  nsweeps = length(time_range) - 1
  szs = zeros(nsweeps, N)
  function region_callback(problem; sweep, kws...)
    return szs[sweep, :] = itn.expect(problem.state, "Sz")
  end

  outputlevel = 0
  inserter_kwargs = (; maxdim=100, cutoff=1E-11)
  updater_kwargs = (; solver=rk4_solver)
  res = ns.applyexp(
    H, psi, time_range; inserter_kwargs, outputlevel, region_callback, updater_kwargs
  )
  #@show norm(res)

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
    for j in 1:N
      szs_ed[sweep, j] =
        scalar(dag(prime(psi, s[j])) * itm.op("Sz", s[j]) * psi) / norm(psi)
    end
  end
  println("\nResult from ED:")
  display(szs_ed)
  #@show norm(psi)

  return nothing
end
