# Instance 2: Lagrangian-relaxation dual

This folder contains the Lagrangian-relaxation instance and its dual-oracle
based solvers.

The APL routine solves the convex minimization form of the dual internally by
minimizing `-d(pi)`, then converts the output back to the original dual
maximization convention.

## Main files

- `run_instance2_lagrangian_experiment.m`: reproducible driver for the legacy
  Instance 2 experiment. It generates the random data, builds the master
  iterates, and runs the selected comparison methods.
- `method_7_apl.m`: accelerated proximal level method applied to the Lagrangian
  dual. It returns the best multiplier, best dual value, reference primal
  value, and run statistics.
- `helper_dual_oracle.m`: evaluates the dual objective and subgradient using Gurobi.
- `helper_primal_reference.m`: computes the primal reference value using Gurobi.
- `method_3_adgpb.m`, `method_4_p_adgpb.m`, and `method_5_pb.m`:
  comparison bundle methods.
- `helper_original_problem.m` and `helper_maxcut_subproblem.m`: auxiliary optimization helpers.

## Paper Mapping

| File | Paper method / role |
| --- | --- |
| `method_7_apl.m` | APL: accelerated proximal-level method for the Lagrangian dual model. |
| `method_3_adgpb.m` | Ad-GPB comparison routine for the dual model. |
| `method_4_p_adgpb.m` | P-Ad-GPB comparison routine initialized with an external stepsize. |
| `method_5_pb.m` | Classical proximal bundle comparison routine. |
| `helper_dual_oracle.m` | Dual objective and subgradient oracle. |
| `helper_primal_reference.m` | Primal reference objective helper. |
| `helper_original_problem.m` | Auxiliary original-problem helper used by legacy bundle code. |
| `helper_maxcut_subproblem.m` | Auxiliary max-cut subproblem helper used by legacy code. |

## Example

Run the full legacy comparison driver:

```matlab
[results, data] = run_instance2_lagrangian_experiment();
```

For a smaller method-selection run after generating the same master iterates:

```matlab
opts = struct( ...
    'methods', {{'APL'}}, ...
    'S_apl', 100, ...
    'apl_opts', struct('max_inner', 40, 'time_limit', 1800));
[results, data] = run_instance2_lagrangian_experiment(opts);
```

The Lagrangian methods can also be called directly: construct or load `q`,
`W`, `T`, `h`, `x`, and an initial multiplier `pi_0`, then compute the primal
reference value and run APL:

```matlab
opt = helper_primal_reference(q, W, T, h, x);
opts = struct('max_inner', 40, 'time_limit', 1800);
[pi, dual_value, primal_ref, stats] = method_7_apl(q, W, T, h, x, pi_0, opt, 100, opts);
```

## Notes

Run MATLAB from this folder, or add this folder to the MATLAB path before
calling the functions. This instance requires a working Gurobi MATLAB
installation. Some helper routines assume the original 120-variable/50-row
experiment dimensions used in the legacy Lagrangian-relaxation tests.
