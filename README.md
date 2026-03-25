# Mohr-Coulomb Material Model for the SIFEL Template

This repository contains a C++ rewrite of the 2D associative, perfectly plastic Mohr-Coulomb constitutive solver from the supplied Matlab code, adapted to the SIFEL-style `vector` and `matrix` framework that was sent in the email thread.

The point of the rewrite was not just to "translate Matlab to C++". The real task was to fit the constitutive algorithm into the way SIFEL expects a material model to work:

- stress is computed in one procedure,
- the stiffness matrix is computed in another procedure,
- intermediate data can be stored in `other`,
- converged history variables are stored in `eqother`,
- `updateval()` transfers converged values at the end of the equilibrium step.

That design comes directly from Tomáš Koudelka's explanation in the emails. The main technical issue raised by Standa was that the Mohr-Coulomb tangent cannot be assembled efficiently from stress alone. The return type, ordered trial eigenvalues, eigenprojections, Hessians, and principal stresses are all needed again in the tangent routine. This repository therefore stores exactly those helper quantities in `other`, so the stiffness routine can reuse them instead of recomputing the whole spectral decomposition.

## What is in the repository

- `matmodel.h`, `matmodel.cpp`
  The material model itself. This is now the main implementation.
- `vector.h`, `vector.cpp`
  The SIFEL-style vector implementation supplied by the author.
- `matrix.h`, `matrix.cpp`
  The matching SIFEL-style matrix implementation.
- `test_matmodel.cpp`
  A regression test written against the framework-level API.
- `matlab_export/`
  The original Matlab reference code, plus the export helper used during validation.

The old standalone kernel split has been removed. There is no separate "core" layer anymore. The constitutive logic now lives directly inside `matmodel.cpp`, because that is the form that is useful for the SIFEL template and closest to the workflow described in the emails.

## Model implemented here

The implementation covers:

- 2D plane strain in 2D-in-3D Voigt ordering
- associative flow rule
- perfect plasticity
- small strains
- elastic response
- return to smooth yield face
- return to left edge
- return to right edge
- return to apex

The component ordering is:

- strain: `[eps_xx, eps_yy, gamma_xy, eps_zz]`
- stress: `[sig_xx, sig_yy, tau_xy, sig_zz]`

The tangent is assembled as a full `4 x 4` matrix for the same ordering.

## How the Matlab code was translated

The Matlab side has two constitutive pieces:

- `constitutive_problem.m`
- `stiffness_matrix.m`

The first routine computes the stress update and all spectral quantities needed by the algorithm. The second one uses those spectral quantities to assemble the consistent tangent.

The C++ rewrite follows the same structure internally, but it is expressed in the idiom expected by the SIFEL template:

1. Read the material parameters into the `matmodel` object.
2. In `nlstresses()`, compute the trial strain from total strain minus plastic strain from `eqother`.
3. Perform the same ordering of trial eigenvalues as in Matlab.
4. Decide the return type.
5. Compute principal stresses, full stress, and updated plastic strain.
6. Pack all helper quantities into the current state buffer `statev`, which plays the role of SIFEL's `other`.
7. In `stiffmat()`, reuse the already prepared `other`-like data instead of rebuilding it from scratch.
8. In `updateval()`, copy the converged history variables into `eqother`.

The important point is that the conversion is not just algebraically equivalent. It follows the storage strategy proposed in the email thread so that the tangent routine has access to the same helper data that was already computed in the stress routine.

## `other` and `eqother`

This is the key design decision in the whole repository.

### `eqother`

`eqother` stores only the quantities that need to survive between equilibrium steps. In this model that means the converged plastic strain:

| index | meaning |
| ---: | --- |
| 0 | `epsp_xx` |
| 1 | `epsp_yy` |
| 2 | `gammap_xy` |
| 3 | `epsp_zz` |

So:

- `ncompeqother = 4`

### `other`

`other` stores the current Newton-step cache:

| range | meaning |
| --- | --- |
| `0..3` | current plastic strain |
| `4` | return type stored as `double` |
| `5..7` | ordered trial principal strains |
| `8..19` | eigenprojections / first derivatives |
| `20..67` | flattened `4 x 4` Hessians in column-major order |
| `68..70` | principal stresses |

So:

- `ncompother = 71`

This layout is exactly the answer to the question raised in the emails: which intermediate quantities need to be passed from the stress procedure to the tangent procedure so that the tangent can be assembled without duplicated work.

## How the current wrapper matches the email workflow

The emails describe the real SIFEL situation: each material point owns its own `other` and `eqother`.

In this standalone repository we do not have the full SIFEL material-point object, but the wrapper follows the same logic:

- `nlstresses()` computes stress and fills the current `statev` buffer
- the same data are cached so that `stiffmat()` can reuse them
- `stiffmat_from_statev()` assembles the tangent directly from the packed current-state buffer
- `updateval()` copies the converged plastic strain from current state to equilibrium state

So the mechanism is the same as the one Tomáš described. In the real SIFEL integration, the cached current-state buffer would live in the integration point's `other` array. Here, the standalone wrapper emulates the same flow inside the template-level class.

## Files worth reading first

If someone needs to understand the code quickly, the best order is:

1. `matmodel.h`
2. `matmodel.cpp`
3. `test_matmodel.cpp`
4. `matlab_export/constitutive_problem.m`
5. `matlab_export/stiffness_matrix.m`

`matmodel.h` shows the public contract and state layout. `matmodel.cpp` shows how the constitutive algorithm was rewritten. The Matlab files make it easy to verify that the algebra and branching logic were preserved.

## Validation

During the earlier replication work, the Matlab model and the C++ rewrite were checked against each other on the final state of the Matlab slope-stability example. The constitutive outputs matched to machine precision. The repository no longer keeps the generated export files, but the helper script remains in `matlab_export/export_loading_process_final.m` if the comparison needs to be repeated.

For day-to-day checking, the repository includes a smaller regression test.

Build:

```bash
make
```

Run:

```bash
./test_matmodel
```

The test covers:

- all five return modes,
- stress computation,
- transfer from current state to converged state,
- consistency of `stiffmat()` with the cached current-state data,
- analytical tangent versus finite differences.

Expected final line:

```text
TEST STATUS: PASS
```

## Notes on the framework files

The `vector.cpp` and `matrix.cpp` sources are treated here as provided framework code. They compile and work for this repository, although the compiler may emit warnings inside those files. The Mohr-Coulomb rewrite was intentionally kept separate from any cleanup of the supplied algebra framework.

## Limitations

This repository implements only the first target discussed in the emails:

- 2D
- associative
- perfectly plastic
- plane strain

It does not implement:

- plane stress,
- hardening or softening,
- a full production SIFEL material-point class with `ipp / im / ido` arguments.

What it does provide is the constitutive logic, the state-variable strategy, and the framework-level wrapper in the exact shape needed to move into the real SIFEL codebase with minimal rewriting.
