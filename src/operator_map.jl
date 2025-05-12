import ITensorNetworks as itn

#
# This function is a workaround for the slow contraction order
# heuristic in ITensorNetworks/src/treetensornetworks/projttns/projttn.jl
# in the projected_operator_tensors(P::ProjTTN) function (line 97 or so)
#

function operator_map(P::itn.ProjTTN, ψ)
  ψ = copy(ψ)
  if itn.on_edge(P)
    for edge in itn.incident_edges(P)
      ψ *= itn.environment(P, edge)
    end
  else
    region = itn.sites(P)
    ie = itn.incident_edges(P)
    # TODO: improvement ideas
    # - check which vertex (first(region) vs. last(region)
    #   has more incident edges and contract those environments first
    # - use automatic contraction order finding
    for edge in ie
      if itn.dst(edge) == first(region)
        #println("Applying E[$edge]")
        ψ *= itn.environment(P, edge)
      end
    end
    for s in itn.sites(P)
      #println("Applying O[$s]")
      ψ *= itn.operator(P)[s]
    end
    for edge in ie
      if itn.dst(edge) != first(region)
        #println("Applying E[$edge]")
        ψ *= itn.environment(P, edge)
      end
    end
  end
  #ITensors.pause()
  return noprime(ψ)
end
