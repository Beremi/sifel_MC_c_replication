# Mohr-Coulomb SIFEL Model

This repository validates SIFEL-shaped Mohr-Coulomb material models in a
standalone C++ harness. The compatibility files here are not a replacement for
SIFEL; they only provide enough of the surrounding interfaces to compile and
test material-point behavior.

Compatibility layer:

- `iotools.h`
- `global.h`
- `mechmat.h`
- `sifel_compat.cpp`

The material-point overloads tested in this repository are:

- `nlstresses(strain, eqstatev, stress, statev)`
- `stiffmat(strain, statev, stress, d)`
- `updateval(statev, eqstatev)`

The SIFEL-style integration-point wrappers are also present and compile against
the local compatibility layer:

- `nlstresses(long ipp, long ido)`
- `matstiff(matrix &d, long ipp, long ido)`
- `updateval(long ipp, long im, long ido)`
- `giveirrstrains(long ipp, long ido, vector &epsp)`

## 2D Implementation

The 2D implementation uses the latest SIFEL-integrated UGN
Mohr-Coulomb files already prepared for SIFEL:

- `mohrc_ugn.h`
- `mohrc_ugn.cpp`

It is a 4-component plane-strain material-point model.

Voigt ordering:

- strain: `[eps_xx, eps_yy, gamma_xy, eps_zz]`
- stress: `[sig_xx, sig_yy, tau_xy, sig_zz]`

Material input read by `read()`:

```text
young poisson cohesion phi
```

`eqstatev`:

| Range | Count | Meaning |
| --- | ---: | --- |
| `0..3` | 4 | converged plastic strain `epsp` |

`statev`:

| Range | Count | Meaning |
| --- | ---: | --- |
| `0..3` | 4 | current plastic strain `epsp` |
| `4..7` | 4 | previous plastic strain `epsp_prev` |
| `8` | 1 | `return_type` stored as `double` |

Total sizes:

- `NCOMP_STRAIN = 4`
- `NCOMP_STRESS = 4`
- `NCOMP_EQOTHER = 4`
- `NCOMP_OTHER = 9`

MATLAB reference:

- `matlab_export/` contains the 2D MATLAB loading-path export.
- `compare_matlab_export.cpp` compares the C++ 2D response against
  `replication_output/matlab/`.

## 3D Implementation

The 3D implementation is added side by side:

- `mohrc3d_ugn.h`
- `mohrc3d_ugn.cpp`

It keeps the same public SIFEL-facing shape as `mohrc_ugn`, but extends the
material point to six strain/stress components and ports the 3D return mapping
from `sysala/slope_stability` commit `8401cc5`
(`+CONSTITUTIVE_PROBLEM/constitutive_problem_3D.m`).

Public SIFEL ordering:

- strain: `[eps_xx, eps_yy, eps_zz, gamma_yz, gamma_xz, gamma_xy]`
- stress: `[sig_xx, sig_yy, sig_zz, tau_yz, tau_xz, tau_xy]`

Internal 3D return-mapping ordering:

```text
[xx, yy, zz, xy, yz, xz]
```

The conversion between public SIFEL ordering and internal/upstream ordering is
done only at the material/API boundary:

- SIFEL to upstream: `[1, 2, 3, 6, 4, 5]`
- upstream to SIFEL: `[1, 2, 3, 5, 6, 4]`

Material input read by `read()` is intentionally the same as 2D:

```text
young poisson cohesion phi
```

`eqstatev`:

| Range | Count | Meaning |
| --- | ---: | --- |
| `0..5` | 6 | converged plastic strain `epsp` in SIFEL order |

`statev`:

| Range | Count | Meaning |
| --- | ---: | --- |
| `0..5` | 6 | current plastic strain `epsp` in SIFEL order |
| `6..11` | 6 | previous plastic strain `epsp_prev` in SIFEL order |
| `12` | 1 | `return_type` stored as `double` |

Total sizes:

- `NCOMP_STRAIN = 6`
- `NCOMP_STRESS = 6`
- `NCOMP_EQOTHER = 6`
- `NCOMP_OTHER = 13`

Return type codes are shared with the 2D model:

| Code | Meaning |
| ---: | --- |
| `0` | elastic |
| `1` | smooth face |
| `2` | left edge, `sigma1 = sigma2` |
| `3` | right edge, `sigma2 = sigma3` |
| `4` | apex, `sigma1 = sigma2 = sigma3` |

MATLAB reference:

- `matlab_export_3d/` intentionally mirrors `matlab_export/` file by file:
  `input_data.m`, `regular_mesh.m`, `preprocessing.m`,
  `constitutive_problem.m`, `stiffness_matrix.m`, `newton.m`,
  `loading_process.m`, and `transformation.m`.
- The 3D differences are contained inside those matching files: tetrahedral
  mesh generation, six strain/stress components, 3D constitutive formulas, and
  the SIFEL/upstream ordering conversions used by export scripts.
- `export_branch_catalog_final.m` writes the five-branch deterministic material
  catalog to `replication_output_3d/`.
- `export_loading_process_final.m` writes the full 3D homogeneous slope export
  to `replication_output_3d/full/`. The default `N_h = 4` dataset covers every
  implemented Mohr-Coulomb return branch.

Additional 3D export documentation and generated plots are in
`matlab_export_3d/README.md`.

## Main 2D/3D Differences

| Area | 2D `mohrc_ugn` | 3D `mohrc3d_ugn` |
| --- | --- | --- |
| Components | 4 plane-strain components | 6 full 3D components |
| Public order | `[xx, yy, xy, zz]` | `[xx, yy, zz, yz, xz, xy]` |
| State size | `eqstatev = 4`, `statev = 9` | `eqstatev = 6`, `statev = 13` |
| MATLAB mesh | 2D triangular slope | 3D P1 tetrahedral slope |
| MATLAB structure | `matlab_export/` | `matlab_export_3d/`, same base file names |
| Reference data | full 2D loading path | branch catalog plus full 3D all-branch path |
| Ordering conversion | none outside 2D layout | SIFEL/upstream conversion at export/API boundary |

## Verification

Build all standalone checks:

```bash
make clean && make
```

Regenerate MATLAB reference exports:

```bash
matlab -batch "cd('matlab_export'); export_loading_process_final; exit"
matlab -batch "cd('matlab_export_3d'); export_branch_catalog_final; export_loading_process_final; exit"
```

Run validation:

```bash
./test_mohrc_ugn
./compare_matlab_export
./test_mohrc3d_ugn
./compare_matlab_export_3d
./compare_matlab_export_3d replication_output_3d
```

Expected results:

```text
MOHRC_UGN TEST STATUS: PASS
COMPARE STATUS: PASS
MOHRC3D_UGN TEST STATUS: PASS
COMPARE 3D STATUS: PASS
```
