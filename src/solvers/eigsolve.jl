using KrylovKit: KrylovKit

function eigsolve_solver(
  operator,
  init;
  which_eigval=:SR,
  ishermitian=true,
  tol=1e-12,
  krylovdim=3,
  maxiter=1,
  verbosity=0,
  eager=false,
  kws...,
)
  howmany = 1
  vals, vecs, info = KrylovKit.eigsolve(
    operator,
    init,
    howmany,
    which_eigval;
    ishermitian,
    tol,
    krylovdim,
    maxiter,
    verbosity,
    eager,
  )
  return vals[1], vecs[1]
end
