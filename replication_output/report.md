# Matlab to C++ Replication Report

## Run Summary

- Final accepted load-step count: 70
- Final printed step index: 71
- Final gravity load factor zeta: 6.82812499999999201e+00
- Final settlement at point A: 3.56949434683096012e+00
- Final Newton iteration count: 29
- Mesh size: nv = 2761, nt = 5300
- Exit reason code: 1 (1 = too_large_settlement, 2 = too_small_load_increment, 3 = max_steps)

## Return Types

| Return type | Matlab count | C++ count |
| --- | ---: | ---: |
| elastic | 4884 | 4884 |
| smooth | 388 | 388 |
| left_edge | 0 | 0 |
| right_edge | 0 | 0 |
| apex | 28 | 28 |

- Return-type mismatches: 0

## Quantity Comparison

| Quantity | Shape | Max abs diff | Mean abs diff | RMS abs diff | Rel L2 diff | Worst location |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| stress S | 4x5300 | 3.183231e-12 | 4.068727e-13 | 6.370813e-13 | 2.676922e-16 | row 2, elem 2701 |
| plastic strain Ep | 4x5300 | 4.440892e-16 | 1.878999e-17 | 2.714993e-17 | 1.569664e-16 | row 1, elem 4998 |
| ordered eig | 3x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| principal stress sigma | 3x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| projection Eig_1 | 4x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| projection Eig_2 | 4x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| projection Eig_3 | 4x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| Hessian EIG_1 reduced | 9x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| Hessian EIG_2 reduced | 9x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| Hessian EIG_3 reduced | 9x5300 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | row 1, elem 1 |
| tangent Sderiv reduced | 9x5300 | 1.164153e-10 | 8.192975e-13 | 6.651568e-12 | 3.093313e-17 | row 4, elem 5161 |

## Exported Data

- Matlab inputs and reference outputs: `replication_output/matlab`
- C++ replay outputs: `replication_output/cpp`
- Full 4x4 C++ Hessians are exported as `hess*_cpp_full.txt`.
- Full 4x4 C++ consistent tangent is exported as `D_cpp_full.txt`.
- Matlab only provides reduced 3x3 Hessians and reduced 3x3 tangent, so the numerical comparison uses those reduced blocks.
