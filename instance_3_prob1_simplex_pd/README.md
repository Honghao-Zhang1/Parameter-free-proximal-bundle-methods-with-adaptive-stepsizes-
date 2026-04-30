# Instance 3: simplex-constrained primal-dual prob1

This folder contains methods for the bounded/simplex-constrained model

```matlab
min ||A*x - b||_1  subject to x >= 0, sum(x) <= D.
```

## Main files

- `run_instance3_sparse_l1_experiment.m`: reproducible sparse benchmark generator and
  comparison driver.
- `method_9_apl.m`: primal-dual accelerated proximal level solver. It returns
  the best point, upper bound, lower bound, and run statistics.
- `helper_project_simplex.m`: projection onto `{x >= 0, sum(x) <= D}`.
- `method_8_adgpb.m`: adaptive proximal bundle comparison method.

## Paper Mapping

The active driver in this folder runs Ad-GPB and the primal-dual APL method
by default. This cleaned folder only includes the Ad-GPB and APL routines for
this instance, not GPB/PB/BG variants.

| File | Paper method / role |
| --- | --- |
| `run_instance3_sparse_l1_experiment.m` | Sparse data generator and experiment driver for the simplex-constrained model. |
| `method_9_apl.m` | APL: accelerated proximal-level method with primal-dual lower-bound updates. |
| `method_8_adgpb.m` | Ad-GPB comparison routine for the same simplex-constrained model. |
| `helper_project_simplex.m` | Helper, not a numbered method; Euclidean projection onto `{x >= 0, sum(x) <= D}`. |

## Example

```matlab
run_instance3_sparse_l1_experiment(1000, 20000, 0.01, 100, 40)
```

To run only one method:

```matlab
opts = struct('methods', {{'APL'}});
run_instance3_sparse_l1_experiment(1000, 20000, 0.01, 100, 40, opts)
```

For a quick smoke test:

```matlab
opts = struct( ...
    'methods', {{'Ad-GPB', 'APL'}}, ...
    'K_bundle', 5, ...
    'S', 2, ...
    'apl_opts', struct('max_inner', 10, 'time_limit', 20));
run_instance3_sparse_l1_experiment(5, 10, 0.2, 2, 40, opts)
```

## Notes

Run MATLAB from this folder, or add this folder to the MATLAB path before
calling the functions. The active APL path is self-contained in this folder.
