# NetworkSolvers.jl

## Package Description

Experimental solvers (`eigsolve` aka `dmrg`, `tdvp`, and others) based on an iterator design.

Design based on DifferentialEquations.jl Integrator Interface:
https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/

Matt's Gist for an early iterators design scheme is here:
https://gist.github.com/mtfishman/fc15f9c675278efb62754b21a1cc7c7e

## To Do List

- Misc. improvements:
    - [ ] replace calls like `tn[v] = T` with `set_index_preservegraph!`
    - [ ] Make root vertex be 1 for MPS / path_graph?
    - [ ] Should 'nsites' keyword be put into something like extracter/updater/inserter kwargs or is 
          it more "global". Can it be left out and just inferred from region sizes?
    - [X] Split subspace arguments like maxdim back into subspace_kwargs, not using truncation_kwargs.
          Rename truncation_kwargs back to inserter_kwargs.
    - [X] Reconceptualize local_tensor as local_state, which might be a collection of tensors, etc.

- Test improvements:
  - [ ] Make examples runnable by runtests.jl
  - [ ] Test quench with 2-site TDVP + subspace

- One-site TDVP redesign:
  - [ ] Fix inserter in 1-site TDVP to apply truncation.
  - [X] Do one-site TDVP in a 2-site "style", just updating the left and bond tensor and not touching
        the 'right' tensor except possibly during the inserter step.
  - [X] How to address issue of "jumps" in sweeping plan?
        How does ITensorNetworks handle it? 
        Just figure out 'next' bond using graph logic?

- [ ] TDVP improvements
  - [ ] Possible bug when combining 1-site and subspace on tree networks
        (came up for a network which was mostly a chain but had an extra ancilla site
         connecting to the impurity site located in middle of the chain)
  - [ ] Implement Euler tour type plan for TDVP. Need to figure out if it 
        does the right amount of time stepping per sweep and what order it is.
  - [X] Fix "densitymatrix" subspace expansion to work with 1-site TDVP and test.
  - [X] Implement 2-site region plan (maybe as separate function for now)
  - [X] Better way of detecting end of sweep and advancing the time step.
        Possibly by querying `isnothing(next_region(region_iterator))`.
  - [X] Change "time" nomenclature to "exponents"? 
        Time can be misleading since there is no "im" included.
        Or put in "im" when calling through `tdvp`, but
        not when calling through `exponentiate` (is `integrate` a better name?)
  - [X] Redesign TDVP around an array of time steps

- DMRG improvements
    - [X] Slow performance for QN-conserving case
    - [X] Alternative sweeping schemes. Go into each subtree then
          back out ("Euler tour"). May help subspace expansion to be more effective.
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


## Timing Notes

- May 11 (main 60cbc1c). Representative timings of DMRG and TDVP.
  ```
  NetworkSolvers DMRG.  nsites=1 
  subspace_kwargs=(; algorithm="densitymatrix", maxdim=4)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 16.1s
  conserve_qns=true:  10.3s

  NetworkSolvers DMRG. nsites=2
  subspace_kwargs=(;)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 10.1s
  conserve_qns=true:  12.1s <-- slower than without QNs

  NetworkSolvers TDVP. (using "test_tdvp.jl")
  total_time=1.0, time_step=0.1, conserve_qns=false, maxdim=16
  nsites=1:  4.2s
  nsites=2: 11.2s

  ITensorMPS DMRG. nsites=2
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 5.1s
  conserve_qns=true:  4.6s
  ```

- After operator map workaround:
  ```
  NetworkSolvers DMRG.  nsites=1 
  subspace_kwargs=(; algorithm="densitymatrix", maxdim=4)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 14.7s
  conserve_qns=true:  10.3s
  --> why are bond dimensions saturating maxdims?
      subspace yes, but is cutoff not being applied during inserter step?

  NetworkSolvers DMRG. nsites=2
  subspace_kwargs=(;)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false:  5.6s
  conserve_qns=true:   5.1s
  --> timings much closer to ITensorMPS DMRG now (bond dimensions slightly different).

  NetworkSolvers TDVP. (using "test_tdvp.jl")
  N=100, total_time=1.0, time_step=0.1, conserve_qns=false, maxdim=16
  nsites=1: 3.2s (down from 4.2s, 1.3x speedup)
  nsites=2: 7.2s (down from 11.2s, 1.5x speedup)
  ```

- After switching to optimal_map (automatic contraction sequence):
  ```
  NetworkSolvers DMRG.  nsites=1 
  subspace_kwargs=(; algorithm="densitymatrix", maxdim=4)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 9.8s
  conserve_qns=true:  7.7s

  NetworkSolvers DMRG. nsites=2
  subspace_kwargs=(;)
  nsweeps=5, N=100, cutoff=1E-9, maxdim=[10, 40, 80, 160]
  conserve_qns=false: 5.8s
  conserve_qns=true:  5.5s 

  NetworkSolvers TDVP. (using "test_tdvp.jl")
  N=100, total_time=1.0, time_step=0.1, conserve_qns=false, maxdim=16
  nsites=1: 3.1s
  nsites=2: 6.4s
  ```


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
