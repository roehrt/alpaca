# Solver Description

## Abstract

We present Alpaca, an exact solver for the Hitting Set and Dominating Set problems, developed for the  [PACE challenge 2025](https://pacechallenge.org/2025).
Both problems are encoded using the straightforward reduction to MaxSAT, followed by instance preprocessing with `maxpre2`.
The resulting formula is then solved using the OLL algorithm, as used in solvers such as `rc2` and `EvalMaxSAT`.

We introduce three enhancements to the standard OLL algorithm: (1) fast start of OLL using a greedy set packing, (2)
determination of zero-cost satisfiability using `kissat`, and (3) a heuristic for unsatisfiable core selection based
on the height of the totalizer tree associated with a literal.

For Hitting Set instances with small set cardinalities, Alpaca uses `SCIP`, benefiting from comparatively stronger LP relaxations.

In our preliminary testing, Alpaca solves all 100 Hitting Set and all 100 Dominating Set public benchmark instances
within the 30 minute time limit per instance.

## MaxSAT Solving

We preprocess all instances using `maxpre2` [^7] with intrinsic at-most-one constraint detection enabled.

Alpaca uses the OLL algorithm for MaxSAT solving [^9] with the Iterative Totalizer encoding [^8] for cardinality constraints
and a simple literal-flipping heuristic for core minimization. These techniques are also used in state-of-the-art MaxSAT solvers
like `rc2` [^6] and `EvalMaxSAT` [^1]. We use `CaDiCaL` [^2] as the underlying SAT solver.

We introduce three improvements over standard OLL:

1. **Set Packing Lower Bound**. A set packing is a valid collection of disjoint albeit trivial unsatisfiable cores for the MaxSAT instance. 
   So before entering the OLL main loop, we run a greedy set packing heuristic and handle these cores. This has two effects:
    - _Fast start_: Although the SAT solver is also able to easily able to identify these trivial cores, a SAT call is costly and is only able to identify one per call.
    - _More balanced totalizer tree_: The SAT solver may find a worse set packing or already tries to reason with the relaxation literals of already relaxed cores.

   Although, after preprocessing, the instance might no longer be a Hitting Set instance, the generalization to MaxSAT is simple: a clause is a trivial
   unsatisfiable core if all its literals are negations of soft literals.

2. **Initial SAT solve**. Problems obtained via reductions to Hitting Set, or to MaxSAT more generally, may exhibit many
   intrinsic at-most-one relationships.
   During presolving, these can simplify the instance into a problem, where the main challenge becomes deciding
   satisfiability or distinguishing between zero-cost and one-cost solutions.
   To accelerate both types of satisfiability checks, we solve the instance (with soft literals hardened) using `kissat` [^3].
   This special-case handling appears to help for many Dominating Set challenge instances.

3.  **Heuristic based on Totalizer Tree height**. MaxSAT solvers usually try to discover multiple, possibly disjoint, unsatisfiable cores
    and minimize them. However, the decision to use a discovered core is often based solely on its size. One important but
    often overlooked factor is the height and size of the tree connected to each literal within the core.
    In Alpaca, we track the height of these trees and prefer cores with literals belonging to shallower trees.
    Besides technical benefits such as improved propagation strength, this heuristic may be intuitively justified by the
    idea that it is easier to reason about literals in shallower trees than those associated with larger subtrees,
    since they usually represent a less abstract property.
    We also apply this heuristic to guide the order in which literals are tried for the literal-flipping heuristic during
    core minimization.
    Interestingly, a [TODO comment](https://github.com/google/or-tools/blob/474b5c337f5a4ed7a7a5e87f1cae2de7dc311b1e/ortools/sat/optimization.cc#L135) in Google's CP-SAT solver indicates that they are aware of this heuristic, at
    least in the context of core minimization.

Additional notes:
- The solver itself is deterministic. However, `maxpre2` is not, since it uses time-based SAT calls.
- For the Hitting Set track, we switch from MaxSAT to integer linear programming using `SCIP` [^4] after thirty seconds if the
  average set cardinality is less than $3.5$. The rationale is that instances with smaller set cardinalities have stronger
  linear relaxations, which ILP solvers like `SCIP` can exploit while MaxSAT solvers cannot.
- In each iteration, we accept all cores whose heuristic values are within a specified tolerance of the best core
  identified. We then apply a greedy set packing heuristic to these accepted cores to ensure that the subset selected
  for processing is disjoint.

## Sketch of correctness

The reduction from Hitting Set to MaxSAT is well known in the literature.
The three modifications introduced to the OLL algorithm do not affect its soundness:
1. The set packing is a valid collection of disjoint unsatisfiable cores.
2. Replacing the SAT solver does not alter the correctness of the algorithm.
3. The selection of which cores to process has no impact on the outcome of OLL.

Our external dependencies, `maxpre2`, `kissat`, `cadical`, and `scip`, are well-established solvers, regularly used and
validated in (Max)SAT and ILP competitions, respectively.

## Potential Enhancements

- **Hyperparameter tuning.** There are many hyperparameters in Alpaca, but they are not fine-tuned (yet) due to lack of time and resources.
- **Portfolio/Special case solvers.** Dedicated solvers have been developed for special cases of the Hitting Set problem,
  such as the 2019 PACE winner WeGotYouCovered for Vertex Cover [^5]. Integrating these solvers could significantly improve performance in those specific cases.
- **Start `SCIP` without `maxpre2` preprocessing.** The MaxSAT preprocessing appears to
  slow down `SCIP`.

## References

[^1]: Avellaneda, Florent. 2020. “A Short Description of the Solver
EvalMaxSAT.” *MaxSAT Evaluation* 8: 364.

[^2]: Biere, Armin, Tobias Faller, Katalin Fazekas, Mathias Fleury, Nils
Froleyks, and Florian Pollitt. 2024a. “CaDiCaL 2.0.” In *Computer Aided
Verification - 36th International Conference, CAV 2024, Montreal, QC,
Canada, July 24-27, 2024, Proceedings, Part I*, edited by Arie Gurfinkel
and Vijay Ganesh, 14681:133–52. Lecture Notes in Computer Science.
Springer.
[https://doi.org/10.1007/978-3-031-65627-9\\7](https://doi.org/10.1007/978-3-031-65627-9\_7).

[^3]: Biere, Armin, Tobias Faller, Katalin Fazekas, Mathias Fleury, Nils
Froleyks, and Florian Pollitt. 2024b. “CaDiCaL, Gimsatul, IsaSAT and Kissat Entering the SAT
Competition 2024.” In *Proc. Of SAT Competition 2024 – Solver, Benchmark
and Proof Checker Descriptions*, edited by Marijn Heule, Markus Iser,
Matti Järvisalo, and Martin Suda, B-2024-1:8–10. Department of Computer
Science Report Series b. University of Helsinki.

[^4]: Bolusani, Suresh, Mathieu Besançon, Ksenia Bestuzheva, Antonia Chmiela,
João Dionísio, Tim Donkiewicz, Jasper van Doornmalen, et al. 2024. “The
SCIP Optimization Suite 9.0.” Technical Report. Optimization Online.
<https://optimization-online.org/2024/02/the-scip-optimization-suite-9-0/>.

[^5]: Hespe, Demian, Sebastian Lamm, Christian Schulz, and Darren Strash. 2020. “WeGotYouCovered: The Winning Solver from the PACE 2019 Challenge,
Vertex Cover Track.” In *Proceedings of the SIAM Workshop on
Combinatorial Scientific Computing, CSC 2020, Seattle, USA, February
11-13, 2020*, edited by H. Martin Bücker, Xiaoye Sherry Li, and
Sivasankaran Rajamanickam, 1–11. SIAM.
<https://doi.org/10.1137/1.9781611976229.1>.

[^6]: Ignatiev, Alexey, António Morgado, and João Marques-Silva. 2019. “RC2:
An Efficient MaxSAT Solver.” *J. Satisf. Boolean Model. Comput.* 11 (1):
53–64. <https://doi.org/10.3233/SAT190116>.

[^7]: Ihalainen, Hannes, Jeremias Berg, and Matti Järvisalo. 2022. “Clause
Redundancy and Preprocessing in Maximum Satisfiability.” In *Proceedings
of the 11th International Joint Conference on Automated Reasoning
(IJCAR ’22)*, 13385:75–94. Lecture Notes in Computer Science. Springer.

[^8]: Martins, Ruben, Saurabh Joshi, Vasco Manquinho, and Inês Lynce. 2014.
“Incremental Cardinality Constraints for MaxSAT.” In *Principles and
Practice of Constraint Programming - 20th International Conference, CP
2014, Lyon, France, September 8-12, 2014. Proceedings*, edited by Barry
O’Sullivan, 8656:531–48. Lecture Notes in Computer Science. Springer.
[https://doi.org/10.1007/978-3-319-10428-7\\39](https://doi.org/10.1007/978-3-319-10428-7\_39).

[^9]: Morgado, António, Carmine Dodaro, and João Marques-Silva. 2014.
“Core-Guided MaxSAT with Soft Cardinality Constraints.” In *Principles
and Practice of Constraint Programming - 20th International Conference,
CP 2014, Lyon, France, September 8-12, 2014. Proceedings*, edited by
Barry O’Sullivan, 8656:564–73. Lecture Notes in Computer Science.
Springer.
[https://doi.org/10.1007/978-3-319-10428-7\\41](https://doi.org/10.1007/978-3-319-10428-7\_41).
