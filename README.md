# Alpaca Hitting Set Solver

This is Alpaca, a student-developed Hitting Set and Dominating Set solver for the [PACE challenge 2025](https://pacechallenge.org/2025).
You can find a short description [here](description.md).

## Build

Build the solver using Docker with
```bash
docker build -t alpaca .
```
and run it using
```bash
docker run -i alpaca < sample.gr
```

## Background

I originally intended Alpaca as a playground just for myself to help me understand the inner workings of an
(unweighted) MaxSAT solver better by writing a toy one myself.
Therefore, the code is
- in a single file.
- super duper messy.
- full of random constants with no deeper meaning
- using lots of naive strategies (e.g. for finding cores and minimizing them)

In short, there are better and more stable tools, and Alpaca was never intended
to be used by others anyway.

## Idea List

... for improving Alpaca as a MaxSAT solver.

- **Minimization by Propagation.** Google's CP-SAT solver uses conflict analysis for minimization in what they call [_minimization by propagation_](https://github.com/google/or-tools/blob/474b5c337f5a4ed7a7a5e87f1cae2de7dc311b1e/ortools/sat/optimization.cc#L58).  
  Their idea is to propagate all literals in the core one after another and, upon conflict, expand the conflict until it consists only of literals from the core. This conflict then becomes the new core.  
  However, this is not possible with typical SAT solver interfaces like the one provided by `CaDiCaL` and thus requires customization.

- **(Conflict-based) Timeouts.** One should replace all fixed loops with loops that use more sophisticated stopping criteria, such as timeouts, the quality of the core according to the heuristic compared to previous cores, or the number of conflicts encountered by the SAT solver (for a more deterministic measure).  
  The goal is to avoid spending an unnecessary amount of time on core discovery or minimization.

- **Restarts.** While already used internally by SAT solvers, MaxSAT solvers may also benefit from restarts. Choosing even a slightly different set of cores can drastically affect the solving time.  
  This also enables another strategy: earlier "fast" passes with more aggressive settings. Often, MaxSAT solvers spend too much time on easy instances by searching for and minimizing many unsatisfiable cores.

- **Local Search.** Short local search phases are expected to have three beneficial effects:
    1. Finding a close upper bound may allow us to perform the final SAT call with hardened assumptions.
    2. It provides a better estimate of the remaining time.
    3. It might benefit the solving process in general due to potential phase-saving effects.
  
  EvalMaxSAT seems to be doing [this](https://github.com/FlorentAvellaneda/EvalMaxSAT/blob/ebd2fbe615858608ee302bc6cbbfc663cb392497/lib/EvalMaxSAT/src/EvalMaxSAT.h#L25).