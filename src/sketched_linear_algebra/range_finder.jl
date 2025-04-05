
"""
  Version of range_finder calling a `sample_from_range` function
  that returns random vectors (range elements) one by one.
"""
function range_finder(
  sample_from_range;
  cutoff=0.0, # experimental feature
  domain_size=typemax(Int),
  max_rank=typemax(Int),
  orthogonal_threshold=1E-12,
  oversample=2,
  range_size=nothing,
  north_pass=2,
)
  range_vectors = []
  (max_rank <= 0) && return range_vectors
  if isnothing(range_size)
    q = sample_from_range()
    qnorm = norm(q)
    (qnorm < orthogonal_threshold) && return range_vectors
    range_vectors = [q / qnorm]
    range_size = prod(size(q))
  end

  max_rank = min(max_rank, range_size, domain_size)
  sketch_size = max_rank + oversample
  sketch_size = min(sketch_size, range_size, domain_size)

  (sketch_size <= 0) && return range_vectors

  for k in (length(range_vectors) + 1):sketch_size
    q = sample_from_range()
    for pass in 1:north_pass
      for qprev in range_vectors
        q = q - dot(qprev, q) * qprev
      end
    end
    qnorm = norm(q)
    (qnorm < orthogonal_threshold) && break
    q /= qnorm
    range_vectors = isempty(range_vectors) ? [q] : push!(range_vectors, q)
    (qnorm < cutoff) && break # experimental cutoff feature
  end

  return range_vectors
end


"""
  Version of range_finder calling a `random_vector` function
  that returns random vectors (domain elements) one by one
  and plugs them into a `linear_map` function.
"""
function range_finder(linear_map, random_vector; domain_size=typemax(Int), kws...)
  v = random_vector()
  vsize = prod(size(v))
  if (domain_size < typemax(Int)) && vsize != domain_size
    error(
      "length of random_vector() output (=$(length(v))) should equal domain_size (=$domain_size)",
    )
  end
  sample_from_range() = linear_map(random_vector())
  return range_finder(sample_from_range; domain_size=vsize, kws...)
end


