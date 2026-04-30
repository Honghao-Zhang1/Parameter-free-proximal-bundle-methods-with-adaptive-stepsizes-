# Parameter-free proximal bundle methods with adaptive stepsizes

This repository contains MATLAB implementations used to test proximal bundle
and accelerated proximal level methods for nonsmooth convex optimization.

The code is organized into three instance folders:

- `instance_1_prob1_nonnegative`: least-absolute-deviation `prob1` models
  with the nonnegative constraint `x >= 0`.
- `instance_2_lagrangian_relaxation`: Lagrangian-relaxation experiments that
  call the dual oracle in `instance_2_lagrangian_relaxation/helper_dual_oracle.m`.
- `instance_3_prob1_simplex_pd`: primal-dual/simplex-constrained `prob1`
  models of the form
  `min ||A*x - b||_1` with nonnegativity/simplex-type constraints.

## Main entry points

- `instance_3_prob1_simplex_pd/run_instance3_sparse_l1_experiment.m` builds a reproducible sparse `prob1` instance and
  runs the primal-dual accelerated proximal level solver.
- `instance_3_prob1_simplex_pd/method_9_apl.m` solves
  `min ||A*x - b||_1` subject to `x >= 0` and `sum(x) <= D`, returning the
  best point, upper bound, lower bound, and run statistics.
- `instance_1_prob1_nonnegative/run_instance1_sparse_l1_experiment.m` builds
  a reproducible sparse nonnegative `prob1` instance and runs the comparison
  methods in that folder.
- `instance_1_prob1_nonnegative/method_7_apl.m` is the accelerated proximal
  level comparison routine for nonnegative `prob1` instances.
- `instance_2_lagrangian_relaxation/method_7_apl.m` applies the accelerated
  proximal level framework to the Lagrangian dual model using
  `instance_2_lagrangian_relaxation/helper_dual_oracle.m`.

The remaining method files are bundle-method variants used for comparison
experiments. Instance 3 only includes its Ad-GPB and APL routines.

## Example

From MATLAB, run a sparse simplex-constrained `prob1` experiment with:

```matlab
cd instance_3_prob1_simplex_pd
run_instance3_sparse_l1_experiment(1000, 20000, 0.01, 100, 40)
```

The arguments are:

- `d`: number of rows in the generated matrix `A`.
- `n`: number of decision variables.
- `s`: sparsity level passed to `sprandn`.
- `D`: radius of the simplex-type feasible set.
- `c`: legacy comparison parameter kept for compatibility with older
  commented experiment blocks.

## Notes

- Some older comparison routines still contain commented legacy calls such as
  `Twocut_subproblem1`; the active code paths use the helpers included here.
- The scripts print progress and timing information to the MATLAB console.
  Newer APL functions also return a `stats` structure with exit information.
- Instance 2 uses Gurobi through MATLAB's `gurobi` interface for its
  mixed-integer oracle and reference problems, including
  `instance_2_lagrangian_relaxation/helper_dual_oracle.m`.
