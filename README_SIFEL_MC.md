# 2D associative perfect-plastic Mohr-Coulomb on the SIFEL vector/matrix framework

## What is in this package

This package keeps the Mohr-Coulomb reimplementation directly on the SIFEL-style `vector` / `matrix` framework supplied by the author. The old standalone C kernel layer has been removed; `matmodel.*` is now the implementation.

Files:

- `vector.h`, `vector.cpp`  
  SIFEL-style vector implementation.
- `matrix.h`, `matrix.cpp`  
  SIFEL-style matrix implementation.
- `matmodel.h`, `matmodel.cpp`  
  Framework-level Mohr-Coulomb material model.
- `test_matmodel.cpp`  
  Regression test using the framework API.
- `matlab_export/`  
  Original Matlab reference code and export helper.
- `Makefile`  
  Builds the framework-level regression test.

## Model scope

- Mohr-Coulomb
- associative flow rule
- perfect plasticity
- small strains
- plane strain in 2D-in-3D ordering
- full 4-component stress response
- full 4x4 consistent tangent

Ordering:

- strain: `[eps_xx, eps_yy, gamma_xy, eps_zz]`
- stress: `[sig_xx, sig_yy, tau_xy, sig_zz]`

## Framework-facing API

The public implementation is the `matmodel` class.

Important methods:

- `read(FILE *in)`
- `print(FILE *out)`
- `nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev)`
- `stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d)`
- `stiffmat_from_statev(const vector &statev, matrix &d)`
- `updateval(const vector &statev, vector &eqstatev)`
- `give_num_of_statev()`
- `give_num_of_eqstatev()`

## State-variable layout

### `eqother`

Only converged plastic strain is kept between equilibrium steps:

| index | meaning |
|---:|---|
| 0 | `epsp_xx` |
| 1 | `epsp_yy` |
| 2 | `gammap_xy` |
| 3 | `epsp_zz` |

So:

- `ncompeqother = 4`

### `other`

The current state buffer stores the plastic strain plus the full spectral cache needed by `stiffmat_from_statev()`:

| range | meaning |
|---|---|
| `0..3` | plastic strain |
| `4` | return type |
| `5..7` | ordered trial principal strains |
| `8..19` | first derivatives / projections |
| `20..67` | flattened 4x4 Hessians, column-major |
| `68..70` | principal stresses |

So:

- `ncompother = 71`

## Relation to the Matlab reference

The Matlab archive still provides the reference constitutive routines:

- `matlab_export/constitutive_problem.m`
- `matlab_export/stiffness_matrix.m`

The C++ implementation in `matmodel.cpp` preserves the same return-map partition:

- elastic
- smooth face
- left edge
- right edge
- apex

## Build and test

Build:

```bash
make
```

Run the framework regression test:

```bash
./test_matmodel
```

The test checks:

- return-type detection,
- stress computation,
- state-variable update,
- analytical tangent against finite differences.

Expected final line:

```text
TEST STATUS: PASS
```

## Matlab-side validation helper

`matlab_export/export_loading_process_final.m` runs the Matlab slope-stability example, captures the final accepted load step, and exports the constitutive data needed for cross-checking against the C++ implementation.

## Important assumptions

1. Associative perfect plasticity only.
2. Friction angle `phi` is in radians.
3. Plane stress is not implemented.
4. The real SIFEL integration can reuse the same `other` / `eqother` layout directly.
