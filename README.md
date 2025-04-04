# NetworkSolvers.jl

## Package Description

Experimental solvers (`eigsolve` aka `dmrg`, `tdvp`, and others) based on an iterator design.

Design based on DifferentialEquations.jl Integrator Interface:
https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/

Matt's Gist for an early iterators design scheme is here:
https://gist.github.com/mtfishman/fc15f9c675278efb62754b21a1cc7c7e

## To Do List

- DMRG improvements
    - [ ] QN subspace expansion (check Andrea's fork of SketchedLinearAlgebra)
    - [X] Subspace expansion support
    - [X] Maxdim etc. as a vector support

- Subspace expansion improvements:
    - [ ] QN subspace support
    - [ ] Subspace kwargs tuple
    - [ ] Better or more automatic handling of expansion size

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
