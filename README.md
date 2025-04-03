# NetworkSolvers.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ITensor.github.io/NetworkSolvers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ITensor.github.io/NetworkSolvers.jl/dev/)
[![Build Status](https://github.com/ITensor/NetworkSolvers.jl/actions/workflows/Tests.yml/badge.svg?branch=main)](https://github.com/ITensor/NetworkSolvers.jl/actions/workflows/Tests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ITensor/NetworkSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ITensor/NetworkSolvers.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Package Description

Experimental solvers (`eigsolve` aka `dmrg`, `tdvp`, and others) based on an iterator design.

Design based on DifferentialEquations.jl Integrator Interface:
https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/

Matt's Gist for an early iterators design scheme is here:
https://gist.github.com/mtfishman/fc15f9c675278efb62754b21a1cc7c7e

## To Do List

- [X] DMRG improvements
    - [X] Subspace expansion support
    - [X] Maxdim etc. as a vector support

- [ ] Add a kind of "checkdone" callback feature. Do we just do it 
    for each algorithm?

- [ ] Keyword arg handling
    - try to "tunnel" as much as possible from user side straight through
      to low levels like inserter, then no complexity can arise in between
    - maybe make "arg packer" functions or @kwdef structs? similar to solvers code

- [ ] SweepIterator design questions

    - [ ] Maybe SweepIterator can just be
        `sweep_iterator(problem, kw_list) = [make_region_iterator(problem, kws) for kws in kw_list]`

    - [ ] Add sweep-level keyword arguments. 
        Where?
        Maybe in Sweeps array, like [(1,kws_sweep1), (2, kws_sweep2), ...]
        - [ ] Think of SweepIterator as just an adapter around region_iterator?
        - [ ] Just an iterator of iterators... also plugging in arguments
            to initialize each sub-iterator

    - [ ] Should SweepIterator actually iterate over whole problem (all sweeps)
        by default?
        Maybe... can think of the current version as adapter like
        `sweep_iterator(repeat_iterator(region_iterator)))`

    - [ ] Is SweepIterator just a fancier version of `Iterators.cycle`?
        I.e. one that lets one detect when each cycle is completed?
        Can we just tell people to use `cycle` if they want this behavior?

- [ ] How to allow `sweep_callback` to update problem
      e.g. to update `current_time` field in TDVPProblem.
      (Have callbacks return problem.)

- [ ] Best way to allow user-defined callbacks? (Alternative to Observer.)

- [ ] TDVP improvements

  - [X] Redesign TDVP around an array of time steps

  - [ ] Change "time" nomenclature to "exponents"? 
    Time can be misleading since there is no "im" included.

  - [ ] Or put in "im" when calling through `tdvp`, but
        not when calling through `exponentiate` (is `integrate` a better name?)

  - [ ] Timers through callbacks. How?
        "Timed callback" adapter. Wraps callback in something like `time = @elapsed f()` 
        and grab time through a closure?

- [ ] Demonstrate iterator adapters, such as "take(iter, n)" that takes
      n steps at each iteration.

## Review of Julia iteration interface

Code such as

```
for item in iter
  # body
end
```

is lowered to

```
next = iterate(iter)
while !isnothing(next)
  (item, state) = next
  # body
  next = iterate(iter, state)
end
```


## Installation instructions

This package resides in the `ITensor/ITensorRegistry` local registry.
In order to install, simply add that registry through your package manager.
This step is only required once.
```julia
julia> using Pkg: Pkg

julia> Pkg.Registry.add(url="https://github.com/ITensor/ITensorRegistry")
```
or:
```julia
julia> Pkg.Registry.add(url="git@github.com:ITensor/ITensorRegistry.git")
```
if you want to use SSH credentials, which can make it so you don't have to enter your Github ursername and password when registering packages.

Then, the package can be added as usual through the package manager:

```julia
julia> Pkg.add("NetworkSolvers")
```

## Examples

````julia
using NetworkSolvers: NetworkSolvers
````

Examples go here.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

