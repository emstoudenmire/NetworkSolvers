using TensorOperations # needed for ITensorNetworks.expect to work properly
using ITensors
import ITensorMPS as itm
import ITensorNetworks as itn
using Printf
using Random
import NetworkSolvers as ns

function tdvp(; N=10, total_time=1.0, time_step=0.1)
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

  tdvp_order = 4
  outputlevel = 0
  nsites = 1
  inserter_kwargs = (; maxdim=16, cutoff=1E-10, normalize=true)
  time_range = 0.0:time_step:total_time
  res_1site = ns.tdvp(H, psi0, time_range; nsites, inserter_kwargs, tdvp_order, outputlevel)

  # Using RK solver
  updater_kwargs = (; solver=ns.runge_kutta_solver, order=4)
  res_rk4 = ns.tdvp(
    H, psi0, 0.0:time_step:total_time; inserter_kwargs, updater_kwargs, outputlevel
  )

  # Using 2-site sweeping scheme
  nsites = 2
  updater_kwargs = (; solver=ns.runge_kutta_solver, order=4)
  res_2site = ns.tdvp(
    H, psi0, time_range; nsites, inserter_kwargs, updater_kwargs, outputlevel
  )

  @show inner(res_1site, res_2site)
  @show inner(res_1site, res_rk4)
  @show inner(res_2site, res_rk4)

  return nothing
end

function test_tdvp(; N=4, total_time=2E-3, time_step=5E-4, tdvp_order=4)
  Random.seed!(1)
  s = itn.siteinds("S=1", N)

  os = itm.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 1/2, "S+", j, "S-", j + 1
    os += 1/2, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, s)
  psi0 = itn.random_mps(s; link_space=8)
  psi0 = itn.truncate(psi0; cutoff=1E-12)

  time_range = 0.0:time_step:total_time
  nsweeps = length(time_range) - 1

  szs_tdvp = zeros(nsweeps, N)
  function sweep_callback(problem; sweep, kws...)
    sz_t = real.(itn.expect(problem.state, "Sz"))
    return szs_tdvp[sweep, :] = sz_t
  end

  outputlevel = 0
  inserter_kwargs = (; maxdim=40, cutoff=1E-10, normalize=true)
  nsites = 2

  #subspace_kwargs = (; algorithm="densitymatrix", maxdim=4)
  subspace_kwargs = (;)

  psi_tdvp = ns.tdvp(
    H,
    copy(psi0),
    time_range;
    nsites,
    inserter_kwargs,
    outputlevel,
    sweep_callback,
    subspace_kwargs,
    tdvp_order
  )
  println("\nResult from TDVP:")
  display(szs_tdvp)
  @show itn.norm(psi_tdvp)


  #
  # Use ED to check
  #
  Hx = prod(H)
  psix = prod(psi0)
  psix /= norm(psix)
  expHx = exp(Hx * (-im*time_step))
  szs_ed = zeros(nsweeps, N)
  for sweep in 1:nsweeps
    psix = noprime(psix * expHx)
    psix /= norm(psix)
    for j in 1:N
      szs_ed[sweep, j] = real(scalar(dag(prime(psix, s[j])) * itm.op("Sz", s[j]) * psix))
    end
  end
  println("\nResult from ED:")
  display(szs_ed)
  @show norm(psix)
  @show abs(scalar(dag(prod(psi_tdvp))*psix))

  println()
  @show norm(szs_ed - szs_tdvp)

  maxerr = 0.0
  err_point = nothing
  for i in 1:size(szs_tdvp, 1), j in 1:size(szs_tdvp, 2)
    err = norm(szs_tdvp[i, j]-szs_ed[i, j])
    if err > maxerr
      maxerr = err
      err_point = (i, j)
    end
  end
  @printf("Largest error (%.3E) at i,j=%d,%d\n", maxerr, err_point[1], err_point[2])
  @printf("   TDVP value = %.10f\n", szs_tdvp[err_point...])
  @printf("     ED value = %.10f\n", szs_ed[err_point...])

  return nothing
end
