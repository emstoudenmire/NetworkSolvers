# NetworkSolvers.jl

## Package Description

Experimental solvers (`eigsolve` aka `dmrg`, `tdvp`, and others) based on an iterator design.

Design based on DifferentialEquations.jl Integrator Interface:
https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/

Matt's Gist for an early iterators design scheme is here:
https://gist.github.com/mtfishman/fc15f9c675278efb62754b21a1cc7c7e

## Notes on New One-Site Design

A challenge of introducing a new vertex to hold the "bond" tensor is:
how do we 'detect' that this is what we did in later code?
In NetworkSolvers, we can keep the old region around.

Ideally we might handle this inside the TreeTensorNetwork or ITensorNetwork
by having a gauge object and using its state to change the behavior of
`orthogonalize` or `gauge_walk` etc.

## To Do List

- Test improvements:
  - [ ] Make examples runnable by runtests.jl
  - [ ] Test quench with 2-site TDVP + subspace

- One-site redesign goals:
  - [ ] When ortho center is on a bond, have a new vertex for the bond matrix.
  - [ ] Don't pass region around anymore. Instead have problem (and possibly state
        inside problem) carry region information like in fitting example.
  - [ ] Check that densitymatrix subspace expansion code works for 1-site TDVP.
  - [ ] Check behavior of bond dimension for 1-site DMRG. 
        How does the subspace expansion interact with 1-site SVD truncation?
  - [ ] Get rid of special edge overload of inserter?
  - [ ] Try other sweeping patterns, like breadth-first versus depth-first.
        Does it make subspace more effective for 1-site?

- [ ] TDVP improvements
  - [ ] Implement 2-site region plan (maybe as separate function for now)
  - [ ] Better way of detecting end of sweep and advancing the time step.
        Possibly by querying `isnothing(next_region(region_iterator))`.
  - [ ] Fix "densitymatrix" subspace expansion to work with 1-site TDVP and test.
  - [ ] Change "time" nomenclature to "exponents"? 
        Time can be misleading since there is no "im" included.
        Or put in "im" when calling through `tdvp`, but
        not when calling through `exponentiate` (is `integrate` a better name?)
  - [ ] Timers through callbacks. How?
        "Timed callback" adapter. Wraps callback in something like `time = @elapsed f()` 
        and grab time through a closure?
  - [X] Redesign TDVP around an array of time steps

- DMRG improvements
    - [X] QN subspace expansion
    - [X] Subspace expansion support
    - [X] Maxdim etc. as a vector support

- [ ] Add a kind of "checkdone" callback feature. Do we just do it 
    for each algorithm?

- Keyword argument handling:
  - [X] How to "mix" keyword argument packs?
        Example: subspace expansion taking `inserter_kwargs`. It works
        but then it causes a dependency between the subspace system and
        inserter system.
        Maybe more "conceptual" kwarg packs, like `truncation_kwargs` ?
  - [ ] How best to "route" orthogonal sets of keyword arguments without
        listing them all?
        Like arguments to alternating_update versus region keyword args:
        does each algorithm need to list all valid arguments? Or just one
        set and assume the other is in the kws... pack? 

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


- [ ] Demonstrate iterator adapters, such as "take(iter, n)" that takes
      n steps at each iteration.

- Subspace expansion improvements:
    - [X] QN subspace support
    - [X] Subspace kwargs tuple
    - [X] Better or more automatic handling of expansion size

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
