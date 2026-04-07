# class-group-tabulation

This repository contains two programs:

- [`clgrp_ell`](src/clgrp_ell_main.c) -
  Uses precomputed fundamental class group data to compute the class group structure of the order of index $\ell$ or $\ell^2$, depending on the Kronecker symbol.

- [`ell_count`](src/ell_count_main.c) -
  Reads the output of `clgrp_ell` alongside precomputed fundamental class group data, and produces a CSV of cumulative counts of various properties.

## Building

Requires Meson, a C/C++ compiler, GMP, PARI (optional), and an MPI implementation.

```bash
meson setup build
meson compile -C build
```

Build options (see [meson.options](meson.options)):

| Option          | Default  | Description                             |
| --------------- | -------- | --------------------------------------- |
| `with_pari`     | `false`  | Enable PARI verification in `clgrp_ell` |
| `sqrtmodp_maxp` | `104729` | Maximum prime for sqrtmodp lookup table |

## Usage

Both programs are run under MPI with the argument order `<D_max> <files> <ell> <folder>`:

```bash
srun -n 256 build/clgrp_ell 137438953472 512 5 ./lmfdb
srun -n 256 build/ell_count 137438953472 512 5 ./lmfdb > ell5.csv
```

Slurm scripts are provided in [scripts](scripts) for primes $3 \leq \ell \leq 31$.

## Results

The results from running all the scripts are in [results](results).
The direct output of `clgrp_ell` is stored locally (may be hosted somewhere accessible in the future).

The folder also contains Typst code to create graphs of the results:

```bash
typst c --root results results/main.typ
```

## Vendored libraries

The class group computation is built on three vendored libraries in [libs](libs), each with local modifications:

- **[liboptarith](https://github.com/maxwellsayles/liboptarith)** -
  Optimised arithmetic operations for 32, 64, and 128-bit integers.
- **[libqform](https://github.com/maxwellsayles/libqform)** -
  Ideal arithmetic in imaginary guadratic number fields.
- **[clgrp](https://github.com/amosunov/clgrp-1.3)** -
  Subroutines utilised for the tabulation of class groups of imaginary quadratic fields.
  Licensed under GPL-2.0-or-later.
