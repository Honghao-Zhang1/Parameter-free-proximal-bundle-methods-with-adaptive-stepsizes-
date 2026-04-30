# Instance 1: nonnegative prob1

This folder contains methods for least-absolute-deviation problems with the
nonnegative feasible region

```matlab
min ||A*x - b||_1  subject to x >= 0.
```

The files here are self-contained MATLAB implementations for this
nonnegative prob1 instance.

## Main files

- `run_instance1_sparse_l1_experiment.m`: main sparse-instance generator
  and experiment driver.
- `method_1_polyak_subgradient.m`: P-Sub*.
- `method_2_gpb.m`: GPB*.
- `method_3_adgpb.m`: Ad-GPB*.
- `method_4_p_adgpb_star.m`: P-Ad-GPB*.
- `method_5_pb.m`: PB*.
- `method_6_bg.m`: BG*.
- `method_7_apl.m`: APL*.
- `helper_project_nonnegative_two_cuts.m`: projection helper used by
  `method_7_apl.m`.

## Paper Mapping

The method order used here follows the routines relevant to this instance
folder:

- P-Sub*
- GPB*
- Ad-GPB*
- P-Ad-GPB*
- PB*
- BG*
- APL*

| File | Paper method / role |
| --- | --- |
| `run_instance1_sparse_l1_experiment.m` | Sparse data generator and experiment driver. |
| `method_1_polyak_subgradient.m` | Pol-Sub: Polyak subgradient baseline. |
| `method_2_gpb.m` | GPB*: fixed-stepsize proximal bundle baseline with a two-cut update. |
| `method_3_adgpb.m` | Ad-GPB*: adaptive proximal bundle method with end-of-cycle stepsize growth. |
| `method_4_p_adgpb_star.m` | P-Ad-GPB*: adaptive proximal bundle method initialized with the Polyak scale. |
| `method_5_pb.m` | PB*: classical proximal bundle method. |
| `method_6_bg.m` | BG*: diagnostic classical bundle benchmark using the generated optimizer. |
| `method_7_apl.m` | APL: accelerated proximal-level comparison baseline. |
| `helper_project_nonnegative_two_cuts.m` | Helper, not a numbered method; projection helper for two affine cuts over `x >= 0`. |

## Example

Run the full experiment with:

```matlab
run_instance1_sparse_l1_experiment(1000, 20000, 0.01)
```

For a quick smoke test that only runs APL:

```matlab
opts = struct( ...
    'methods', {{'apl'}}, ...
    'S_apl', 2, ...
    'apl_opts', struct('max_inner', 10, 'time_limit', 20));
run_instance1_sparse_l1_experiment(5, 10, 0.2, opts)
```

## Notes

Run MATLAB from this folder, or add this folder to the MATLAB path before
calling the functions. The legacy methods print their own progress and timing
information; `run_instance1_sparse_l1_experiment.m` also returns a `results` structure for the
baseline and APL* run.
