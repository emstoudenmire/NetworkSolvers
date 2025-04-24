import ITensors as it
import ITensorNetworks as itn
import NamedGraphs as ng
using Graphs: Graphs

default_maxdim() = typemax(Int)
default_mindim() = 1
default_cutoff() = 0.0

get_or_last(x, i::Integer) = (i >= length(x)) ? last(x) : x[i]
