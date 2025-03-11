using Literate: Literate
using NetworkSolvers: NetworkSolvers

Literate.markdown(
  joinpath(pkgdir(NetworkSolvers), "examples", "README.jl"),
  joinpath(pkgdir(NetworkSolvers), "docs", "src");
  flavor=Literate.DocumenterFlavor(),
  name="index",
)
