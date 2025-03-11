using Literate: Literate
using NetworkSolvers: NetworkSolvers

Literate.markdown(
  joinpath(pkgdir(NetworkSolvers), "examples", "README.jl"),
  joinpath(pkgdir(NetworkSolvers));
  flavor=Literate.CommonMarkFlavor(),
  name="README",
)
