using TensorOperations # needed for ITensorNetworks.expect to work properly
using ITensors
import ITensorMPS as mps
using ITensorMPS: OpSum, op
import ITensorNetworks as itn
using Printf
using Random
import NetworkSolvers as ns

include("../test/utilities/simple_ed_methods.jl")
include("../test/utilities/write_data.jl")

function print_linkdims(psi)
  ld = itn.linkdims(psi)
  for e in itn.edges(psi)
    println("  $e: ", ld[e])
  end
end

function quench(;
  N=8, total_time=4.0, time_step=0.05, cutoff=1E-14, maxdim=5000, tdvp_order=4
)
  sites = itn.siteinds("S=1/2", N)

  os = mps.OpSum()
  for j in 1:(N - 1)
    os += "Sz", j, "Sz", j + 1
    os += 0.5, "S+", j, "S-", j + 1
    os += 0.5, "S-", j, "S+", j + 1
  end
  H = itn.mpo(os, sites)

  g = itn.underlying_graph(sites)
  state = Dict{Int,String}()
  for (j, v) in enumerate(itn.vertices(sites))
    state[v] = isodd(j) ? "Up" : "Dn"
  end
  psi0 = itn.mps(state, sites)

  time_range = 0.0:time_step:total_time

  sz0 = real(only(itn.expect(psi0, "Sz", N÷2)))
  szs = [sz0]
  function sweep_callback(problem; kws...)
    sz = real(only(itn.expect(problem.state, "Sz", N÷2)))
    push!(szs, sz)
  end

  outputlevel = 0
  nsites = 1
  truncation_kwargs = (; maxdim, cutoff, normalize=true)
  nsites = 2
  subspace_kwargs = (; algorithm="densitymatrix", expansion_factor=1.2, max_expand=4)
  updater_kwargs = (; solver=ns.runge_kutta_solver, order=4)

  println("Calling TDVP")
  psif = ns.tdvp(
    H,
    psi0,
    time_range;
    nsites,
    sweep_callback,
    subspace_kwargs,
    truncation_kwargs,
    updater_kwargs,
    outputlevel,
    tdvp_order,
  )

  print_linkdims(psif)

  fname = "szs_tdvp_N$N.dat"
  println("Writing file \"$fname\"")
  write_data(fname, time_range, szs)

  if N <= 8
    println("Using ED to check")
    psix = ed_time_evolution(H, psi0, time_range; normalize=true)
    fidelity = abs(scalar(dag(psix)*prod(psif)))
    @printf("Fidelity <psi_exact|psi_tdvp> = %.12f\n", fidelity)
  end

  return nothing
end

function tebd_quench(; N=8, cutoff=1E-8, total_time=4.0, time_step=0.05)
  s = mps.siteinds("S=1/2", N; conserve_qns=true)

  gates = ITensor[]
  for j in 1:(N - 1)
    s1 = s[j]
    s2 = s[j + 1]
    hj =
      op("Sz", s1) * op("Sz", s2) +
      1 / 2 * op("S+", s1) * op("S-", s2) +
      1 / 2 * op("S-", s1) * op("S+", s2)
    Gj = exp(-im * time_step / 2 * hj)
    push!(gates, Gj)
  end
  append!(gates, reverse(gates))

  psi = mps.MPS(s, n -> isodd(n) ? "Up" : "Dn")

  c = N÷2
  szs = Float64[]
  time_range = 0.0:time_step:total_time
  for t in time_range
    Sz = real(only(mps.expect(psi, "Sz"; sites=c)))
    push!(szs, Sz)
    t≈total_time && break
    psi = mps.apply(gates, psi; cutoff)
    normalize!(psi)
  end

  fname = "szs_tebd_N$N.dat"
  println("Writing file \"$fname\"")
  write_data(fname, time_range, szs)
end
