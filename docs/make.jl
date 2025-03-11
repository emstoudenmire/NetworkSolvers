using NetworkSolvers: NetworkSolvers
using Documenter: Documenter, DocMeta, deploydocs, makedocs

DocMeta.setdocmeta!(
  NetworkSolvers, :DocTestSetup, :(using NetworkSolvers); recursive=true
)

include("make_index.jl")

makedocs(;
  modules=[NetworkSolvers],
  authors="ITensor developers <support@itensor.org> and contributors",
  sitename="NetworkSolvers.jl",
  format=Documenter.HTML(;
    canonical="https://ITensor.github.io/NetworkSolvers.jl",
    edit_link="main",
    assets=String[],
  ),
  pages=["Home" => "index.md", "Reference" => "reference.md"],
)

deploydocs(;
  repo="github.com/ITensor/NetworkSolvers.jl", devbranch="main", push_preview=true
)
