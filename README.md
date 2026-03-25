# Mohr-Coulomb Rewrite for the SIFEL API

This repository contains a C++ rewrite of the 2D associative, perfectly plastic Mohr-Coulomb constitutive model from the Matlab reference code. The target form is the SIFEL material-model API based on `vector` and `matrix`, with the constitutive update split into:

- stress computation `S`
- consistent tangent computation `DS`

The implementation is plane strain, small strain, associative, and perfectly plastic. It covers elastic response, return to the smooth yield surface, left edge, right edge, and apex.

## Interface

The public entry points are in [matmodel.h](/home/beremi/repos/sifel_MC_c_replication/matmodel.h) and implemented in [matmodel.cpp](/home/beremi/repos/sifel_MC_c_replication/matmodel.cpp):

- `nlstresses(strain, eqstatev, stress, statev)`
  Computes the current stress vector `S` and fills the current-state buffer `statev`.
- `stiffmat(strain, eqstatev, stress, d)`
  Computes the consistent tangent matrix `DS`.
- `updateval(statev, eqstatev)`
  Transfers converged history variables to the equilibrium-state buffer.

Voigt ordering:

- strain: `[eps_xx, eps_yy, gamma_xy, eps_zz]`
- stress: `[sig_xx, sig_yy, tau_xy, sig_zz]`

## Implementation Notes

The Matlab constitutive logic is split here the same way it is used by the SIFEL API:

- `compute_response()` performs the return mapping and prepares all auxiliary quantities needed later by the tangent routine.
- `compute_tangent()` assembles the consistent tangent from the data already prepared by the stress update.

The C++ rewrite uses the supplied SIFEL algebra layer directly for vector and matrix work. Generic operations are expressed with functions such as `copyv`, `rcopyv`, `nullv`, `addmultv`, `copym`, `nullm`, `addm`, `addmultm`, and `vxv`.

## State Data Passed From `S` to `DS`

The key part of this rewrite is the layout of `statev` (`other`) and `eqstatev` (`eqother`).

`eqstatev` stores only the converged history variables needed to start the next return mapping. In this model that is just the plastic strain from the last equilibrium step.

`statev` stores:

- the updated plastic strain for the current Newton iterate
- the detected return type
- the ordered trial principal strains
- the eigenprojections
- the Hessians of the ordered eigenvalues
- the principal stresses

This avoids recomputing the spectral decomposition and return classification when `stiffmat()` is called after `nlstresses()` for the same material point state.

### `eqstatev` layout

| Index | Symbol | Meaning |
| ---: | --- | --- |
| 0 | `epsp_xx` | plastic strain in `xx` |
| 1 | `epsp_yy` | plastic strain in `yy` |
| 2 | `gammap_xy` | plastic engineering shear strain |
| 3 | `epsp_zz` | plastic strain in `zz` |

Size: `4`

### `statev` layout

| Index range | Symbol | Meaning | Used by `DS` |
| --- | --- | --- | --- |
| `0..3` | `epsp` | current plastic strain in Voigt order | no |
| `4` | `return_type` | return mode stored as `double` | yes |
| `5..7` | `eig(1:3)` | ordered trial principal strains | yes |
| `8..11` | `proj_1` | first eigenprojection / first derivative of `eig_1` | yes |
| `12..15` | `proj_2` | first eigenprojection / first derivative of `eig_2` | yes |
| `16..19` | `proj_3` | first eigenprojection / first derivative of `eig_3` | yes |
| `20..35` | `hess_1` | `4 x 4` Hessian of `eig_1`, flattened | yes |
| `36..51` | `hess_2` | `4 x 4` Hessian of `eig_2`, flattened | yes |
| `52..67` | `hess_3` | `4 x 4` Hessian of `eig_3`, flattened | yes |
| `68..70` | `sigma(1:3)` | principal stresses after return mapping | yes |

Size: `71`

### Hessian storage

Each Hessian is stored as a flattened `4 x 4` block. The packing order is the same one used in [matmodel.cpp](/home/beremi/repos/sifel_MC_c_replication/matmodel.cpp): the inner loop runs over row index `i`, the outer loop over column index `j`.

That means the stored sequence is:

`H(0,0), H(1,0), H(2,0), H(3,0), H(0,1), ..., H(3,3)`

`stiffmat()` unpacks this layout back into a `4 x 4` matrix before tangent assembly.

## Control Flow

For one material-point evaluation, the intended sequence is:

1. `nlstresses()` computes `S` and fills `statev`.
2. `stiffmat()` reuses the cached `statev` content to build `DS`.
3. `updateval()` copies the converged plastic strain from `statev` to `eqstatev`.

If `stiffmat()` is called immediately after `nlstresses()` with the same inputs, the expensive return-mapping part is not recomputed.

## Files

- [matmodel.h](/home/beremi/repos/sifel_MC_c_replication/matmodel.h): public API, constants, state layout
- [matmodel.cpp](/home/beremi/repos/sifel_MC_c_replication/matmodel.cpp): constitutive update and tangent assembly
- [test_matmodel.cpp](/home/beremi/repos/sifel_MC_c_replication/test_matmodel.cpp): regression test
- [matlab_export/](/home/beremi/repos/sifel_MC_c_replication/matlab_export): original Matlab reference implementation
- [vector.h](/home/beremi/repos/sifel_MC_c_replication/vector.h), [vector.cpp](/home/beremi/repos/sifel_MC_c_replication/vector.cpp): supplied vector infrastructure
- [matrix.h](/home/beremi/repos/sifel_MC_c_replication/matrix.h), [matrix.cpp](/home/beremi/repos/sifel_MC_c_replication/matrix.cpp): supplied matrix infrastructure

## Build and Test

```bash
make
./test_matmodel
```

Expected result:

```text
TEST STATUS: PASS
```
