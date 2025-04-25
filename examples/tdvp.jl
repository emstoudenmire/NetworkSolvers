using TensorOperations # needed for ITensorNetworks.expect to work properly
using ITensors
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
using Random
import NetworkSolvers as ns

function tdvp(; N=10, total_time=-1.0, time_step=-0.1)
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

  outputlevel = 0
  nsites = 1
  truncation_kwargs = (; maxdim=16, cutoff=1E-10, normalize=true)
  time_range = 0.0:time_step:total_time
  res_1site = ns.tdvp(H, psi0, time_range; nsites, truncation_kwargs, outputlevel)

  # Using RK solver
  updater_kwargs = (; solver=ns.runge_kutta_solver, order=4)
  res_rk4 = ns.tdvp(
    H, psi0, 0.0:time_step:total_time; truncation_kwargs, updater_kwargs, outputlevel
  )

  # Using 2-site sweeping scheme
  nsites = 2
  res_2site = ns.tdvp(H, psi0, time_range; nsites, truncation_kwargs, outputlevel)

  @show inner(res_1site, res_2site)
  @show inner(res_1site, res_rk4)
  @show inner(res_2site, res_rk4)

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
  truncation_kwargs = (; maxdim=100, cutoff=1E-10)
  updater_kwargs = (; solver=ns.runge_kutta_solver, order=2)
  subspace_kwargs = (; algorithm="densitymatrix", maxdim=4)
  #subspace_kwargs = (;)
  res = ns.applyexp(
    H,
    psi,
    time_range;
    truncation_kwargs,
    outputlevel,
    region_callback,
    subspace_kwargs,
    updater_kwargs,
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

  println()
  @show norm(szs_ed - szs)

  maxerr = 0.0
  err_point = nothing
  for i in 1:size(szs, 1), j in 1:size(szs, 2)
    err = norm(szs[i, j]-szs_ed[i, j])
    if err > maxerr
      maxerr = err
      err_point = (i, j)
    end
  end
  @printf("Largest error (%.3E) at i,j=%d,%d\n", maxerr, err_point[1], err_point[2])
  @printf("   TDVP value = %.10f\n", szs[err_point...])
  @printf("     ED value = %.10f\n", szs_ed[err_point...])

  return nothing
end
