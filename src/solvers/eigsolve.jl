using KrylovKit: KrylovKit

function eigsolve_solver(
  operator,
  init;
  which_eigval=:SR,
  ishermitian=true,
  tol=1e-10,
  krylovdim=40,
  maxiter=2,
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
