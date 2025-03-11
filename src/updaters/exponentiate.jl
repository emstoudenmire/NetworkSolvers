using KrylovKit: KrylovKit

function exponentiate_updater(
  operator,
  init;
  time_step,
  krylovdim=30,
  maxiter=100,
  verbosity=0,
  tol=1E-12,
  ishermitian=true,
  issymmetric=true,
  eager=true,
)
  result, exp_info = KrylovKit.exponentiate(
    operator,
    time_step,
    init;
    eager,
    krylovdim,
    maxiter,
    verbosity,
    tol,
    ishermitian,
    issymmetric,
  )
  return result
end
