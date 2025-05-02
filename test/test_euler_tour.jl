import Graphs: dst, edges, src, vertices
import NetworkSolvers as ns
using Test: @test, @testset

include("utilities/tree_graphs.jl")

@testset "Euler tour" begin
  g = build_tree(; nbranch=3, nbranch_sites=3)

  start_vertex = (0,0)

  edge_tour = ns.euler_tour_edges(g, start_vertex)
  for j in 1:length(edge_tour)-1
    @test dst(edge_tour[j]) == src(edge_tour[j+1])
  end
  for e in edges(g)
    @test e in edge_tour
    @test reverse(e) in edge_tour
  end

  vertex_tour = ns.euler_tour_vertices(g, start_vertex)
  for v in vertices(g)
    #println(v," ",count(==(v),vertex_tour))
    @test v in vertex_tour
  end
end
