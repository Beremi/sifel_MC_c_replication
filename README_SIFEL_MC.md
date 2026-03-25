# 2D associative perfect-plastic Mohr-Coulomb for SIFEL

## What is in this package

This package contains a clean rewrite of the Matlab constitutive kernel from the supplied archive into a small C-style computational core plus a thin wrapper matching the `matmodel.h/.cpp` template that was sent in the email thread.

Files:

- `mcppp2d_core.h`, `mcppp2d_core.cpp`  
  Standalone constitutive kernel. No SIFEL dependency.
- `matmodel.h`, `matmodel.cpp`  
  Template-level wrapper using SIFEL `vector` / `matrix` style.
- `test_mcppp2d.cpp`  
  Standalone regression test.
- `Makefile`  
  Builds the standalone regression test.

## What the code implements

Model scope:

- Mohr-Coulomb
- associative flow rule
- perfect plasticity (no hardening / softening)
- small strains
- closed-form stress return as in the Matlab code
- full 4-component constitutive response in the Voigt ordering
  - strain: `[eps_xx, eps_yy, gamma_xy, eps_zz]`
  - stress: `[sig_xx, sig_yy, tau_xy, sig_zz]`
- full **4x4** consistent tangent matrix for the same ordering

This is the natural constitutive setting for **plane strain** in SIFEL.

## Relation to the Matlab files

The Matlab archive had two essential constitutive pieces:

1. `constitutive_problem.m`  
   Computes stress, ordered trial principal strains, eigenprojections, second derivatives and principal stresses.
2. `stiffness_matrix.m`  
   Computes the consistent tangent from the already known spectral data.

The core rewrite preserves the same return-type logic:

- elastic
- return to smooth yield face
- return to left edge
- return to right edge
- return to apex

The main improvement in this rewrite is that the tangent is extended from the Matlab `3x3` in-plane form to the full `4x4` form needed by the SIFEL template ordering.

## SIFEL state-variable layout

The email discussion makes the `other` / `eqother` split the key design point.

### `eqother` (converged state, used as input to next load step)

`eqother` contains only the variables that must survive the equilibrium-step transition:

| index | meaning |
|---:|---|
| 0 | plastic strain `epsp_xx` |
| 1 | plastic strain `epsp_yy` |
| 2 | plastic strain `gammap_xy` |
| 3 | plastic strain `epsp_zz` |

So:

- `ncompeqother = 4`

### `other` (current Newton-iteration cache)

`other` contains the current plastic strain **plus** every auxiliary quantity needed by the tangent operator, so `matstiff()` can avoid recomputing the entire spectral decomposition.

| range | meaning |
|---|---|
| `0..3` | current plastic strain `[epsp_xx, epsp_yy, gammap_xy, epsp_zz]` |
| `4` | return type as `double` (`0 elastic, 1 smooth, 2 left edge, 3 right edge, 4 apex`) |
| `5..7` | ordered trial principal strains `[eig_1, eig_2, eig_3]` |
| `8..11` | projection `Eig_1` (4 components) |
| `12..15` | projection `Eig_2` (4 components) |
| `16..19` | projection `Eig_3` (4 components) |
| `20..35` | Hessian `EIG_1` as flattened `4x4` matrix |
| `36..51` | Hessian `EIG_2` as flattened `4x4` matrix |
| `52..67` | Hessian `EIG_3` as flattened `4x4` matrix |
| `68..70` | principal stresses `[sigma_1, sigma_2, sigma_3]` |

So:

- `ncompother = 71`

### Hessian flattening convention

The Hessians are packed in **column-major** order to stay close to the original Matlab vectorization:

```text
flat_index = row + 4*col
```

for `row, col = 0..3`.

## Why this split is the right one for SIFEL

Only plastic strain is a true history variable for the next equilibrium step in this perfect-plastic model. Everything else is a current-iteration cache. Therefore:

- `eqother` should stay small and contain only the converged plastic strain,
- `other` should contain the full tangent cache.

That matches the intent from the email thread exactly: `nlstresses()` fills the cache, `matstiff()` reuses it, and `updateval()` copies only the converged history variables.

## Files and API details

## `mcppp2d_core.*`

The computational kernel exposes these functions:

- `mcppp2d_init_params(...)`
- `mcppp2d_validate_params(...)`
- `mcppp2d_compute_response(...)`
- `mcppp2d_compute_tangent(...)`
- `mcppp2d_pack_other(...)`
- `mcppp2d_unpack_other(...)`
- `mcppp2d_pack_eqother(...)`

That is the recommended layer to reuse inside the actual SIFEL material class.

## `matmodel.*`

This wrapper is targeted at the exact template from the email. It adds:

- parameter reading/printing,
- `give_num_of_statev()` and `give_num_of_eqstatev()`,
- `updateval()`,
- `stiffmat_from_statev()`.

Important note: the supplied template signature for

```cpp
stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d)
```

does **not** pass the current `statev` / `other` buffer into the function. Therefore the template wrapper recomputes the cache internally.

For the **real SIFEL material class** this should be changed back to the intended workflow:

- `nlstresses()` stores the cache to `Mm->ip[ipp].other[ido + ...]`
- `matstiff()` reads that cache from `other`
- `updateval()` copies only plastic strain to `eqother`

## How to integrate into the actual SIFEL repository

The exact registration points depend on the checkout, but the implementation pattern should be this:

### 1. Add the new files to `SIFEL/MEFEL/SRC`

Add at least:

- `mcppp2d_core.h`
- `mcppp2d_core.cpp`
- a new SIFEL material-model wrapper, e.g. `mcppp2d.h`, `mcppp2d.cpp`

I recommend **not** using the generic filename `matmodel.*` in the real repository. Use a distinct class/file name to avoid collisions.

### 2. Include the new header in the MEFEL umbrella header list

Add the new material header next to the other material-model headers (the public Doxygen index shows that SIFEL aggregates them in `seqfilesm.h`).

### 3. Register the material type in the material factory / dispatch

Follow the same pattern as existing plastic models such as `j2flow2` and `mohrcoulombpar`.

You will need dispatcher branches for at least:

- reading material data,
- `nlstresses(...)`,
- `matstiff(...)`,
- `updateval(...)`,
- `givencompother(...)`,
- `givencompeqother(...)`.

### 4. Implement the SIFEL wrapper methods

The real SIFEL wrapper should have methods in the same spirit as the existing models:

```cpp
void mcppp2d::nlstresses(long ipp, long im, long ido)
void mcppp2d::matstiff(matrix &d, long ipp, long ido)
void mcppp2d::updateval(long ipp, long im, long ido)
```

### 5. `nlstresses()` logic

Pseudo-code:

```cpp
void mcppp2d::nlstresses(long ipp, long im, long ido)
{
    double strain[4];
    double eqother[4] = {0.0, 0.0, 0.0, 0.0};
    double other[71];
    mcppp2d_cache cache;

    // read current strain from Mm->ip[ipp].strain[0..3]
    // read previous plastic strain from Mm->ip[ipp].eqother[ido + 0..3]

    mcppp2d_compute_response(&par, strain, eqother, &cache);
    mcppp2d_pack_other(&cache, other);

    // write stress to Mm->ip[ipp].stress[0..3]
    // write other[0..70] to Mm->ip[ipp].other[ido + ...]
}
```

### 6. `matstiff()` logic

Pseudo-code:

```cpp
void mcppp2d::matstiff(matrix &d, long ipp, long ido)
{
    mcppp2d_cache cache;
    double other[71];
    double dd[4][4];

    if (Mp->nlman->stmat == initial_stiff)
    {
        Mm->elmatstiff(d, ipp);
        return;
    }

    // read other[0..70] from Mm->ip[ipp].other[ido + ...]
    mcppp2d_unpack_other(other, &cache);
    mcppp2d_compute_tangent(&par, &cache, dd);

    // copy dd to matrix d
}
```

### 7. `updateval()` logic

Only copy the plastic strain to `eqother`:

```cpp
void mcppp2d::updateval(long ipp, long im, long ido)
{
    for (long i=0; i<4; i++)
        Mm->ip[ipp].eqother[ido+i] = Mm->ip[ipp].other[ido+i];
}
```

Do **not** copy the whole `other` block unless you deliberately want the larger memory footprint.

### 8. Tell SIFEL how many buffer components the model needs

Use:

```cpp
ncompother   = 71;
ncompeqother = 4;
```

### 9. Rebuild SIFEL

The public install page states that SIFEL uses the GNU make build system. Add the new `.cpp` file to the appropriate makefile / project file in the checkout and rebuild the mechanical module in the same way as the existing material models.

## Standalone build and test

Build the standalone regression test:

```bash
make
./test_mcppp2d
```

or directly:

```bash
g++ -O2 -std=c++11 -Wall -Wextra -pedantic \
    mcppp2d_core.cpp test_mcppp2d.cpp -o test_mcppp2d
./test_mcppp2d
```

Expected outcome:

- all five return modes are detected correctly,
- the analytical tangent matches a central finite-difference tangent,
- the `other` / `eqother` packing roundtrips correctly,
- final line is `TEST STATUS: PASS`.

## Recommended SIFEL functionality test

The cleanest FE-level verification is a **uniform-strain patch test**.

### Geometry

Use a unit square in plane strain, e.g. one quadrilateral or two triangles.

### Prescribed displacement field

For a target strain state

```text
[eps_xx, eps_yy, gamma_xy, eps_zz]
```

prescribe on the unit square the displacement field

```text
u_x(x,y) = eps_xx*x + 0.5*gamma_xy*y
u_y(x,y) = eps_yy*y + 0.5*gamma_xy*x
```

This generates a constant in-plane strain

```text
eps_xx, eps_yy, gamma_xy
```

throughout the element.

### Suggested loading cases

Use zero initial plastic strain and test these cases in separate runs:

1. Elastic

```text
eps = [ 1.0e-5,  0.0,    0.0,   0.0 ]
```

2. Smooth-face return

```text
eps = [-5.0e-3,  0.0,    0.0,   5.0e-3]
```

3. Left-edge return

```text
eps = [-5.0e-3,  5.0e-3, 0.0,   5.0e-3]
```

4. Right-edge return

```text
eps = [-5.0e-3, -5.0e-3, 0.0,   1.0e-2]
```

5. Apex return

```text
eps = [ 0.0,     0.0,    0.0,   5.0e-3]
```

### What to compare

For each case, compare at the integration point:

- `stress[0..3]`,
- return type stored in `other[4]`,
- plastic strain stored in `other[0..3]`,
- consistent tangent from `matstiff()`.

The easiest reference is the standalone executable `test_mcppp2d`, which prints the stress and checks the tangent numerically.

## Important limitations / assumptions

1. This implementation is for the **associative**, **perfectly plastic** model only.
2. It is intended for the `4`-component strain/stress ordering from the supplied SIFEL template.
3. It does **not** implement plane-stress iteration.
4. It assumes friction angle `phi` is in **radians**.
5. The computational core and the standalone test were compiled and executed here.
6. The repository-integration instructions were written against the supplied template plus the publicly indexed SIFEL interface pages, but I could not fetch the full public Doxygen/SVN pages during this session. So the exact dispatcher filenames in your checkout should be verified locally before the final merge.

## Bottom line

The important architectural decision is:

- keep **plastic strain only** in `eqother`,
- keep the **full tangent cache** in `other`.

That gives the exact workflow requested in the email thread and avoids duplicated spectral work in the real SIFEL `matstiff()` implementation.
