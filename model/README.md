# model/ - the Bubble Model (C)

Simulates the lesions of one macaque lung from infection to day 80, N times (N independent
simulated animals), and writes the daily mean and SD of the CT outputs over the N simulations.
Simulations are distributed over MPI processes and OpenMP threads.

## Build

```sh
make            # macOS: mpicc + Homebrew libomp  (Linux: make OMP=-fopenmp)
make hpc        # Intel oneAPI (mpiicx), as used for the paper's runs
```

Needs an MPI C compiler with OpenMP and C23 support. Lesion size thresholds can be changed at
compile time, e.g. `make DEFS='-DCONS_MIN_MAINAXIS=5.5 -DMICRONOD_MIN_MAINAXIS=0.5'`
(defaults 4.5 and 1 mm; used by the `thresholds` run).

## Run

The model reads and writes relative paths, so it runs inside a directory with this layout
(`runs/launch.py` creates it for every run):

```
config/params_bm.txt              model parameters (below)
config/experimental_results.txt   targets for the error (data/processed/experimental_results.txt)
config/macaque/                   the virtual lung (link to lung/macaque/, with D.bin unpacked)
results/                          outputs
```

```sh
mpirun -n 1 /path/to/model/bubble-model -noplot -dlobe -N 10000 -seed 1
```

- `-N <n>`: number of simulated animals (default 100; the paper uses 10^4).
- `-dlobe`: lesion radius limited by the distance to the lobe wall (used in the paper); without it,
  by the segment wall (`dseg` case of the experiments).
- `-seed <s>`: seed of the random numbers; the same seed gives the same results, for any number of
  processes and threads.
- `-noplot`: no plot window at the end.
- `-server <dir>`: stay running and answer evaluation requests in `<dir>` (used by `runs/launch.py`;
  avoids reading the lung for every evaluation).

The remaining options in `src/main.c` (`-debug_seq`, `-debug_par`, `-pin`, `-video`,
`-write_sample`, `-append`, `-human`) are development options not used in the paper (`-human`
needs a human lung, not included).

## Parameters (`config/params_bm.txt`)

19 lines `name = value`, in this order (read with fixed formats):

| name | meaning | paper |
|---|---|---|
| `v0` | lesion growth rate (day^-1) | fitted |
| `fvimm` | factor on the growth rate after t_imm | fitted |
| `r0` | mean maximum lesion radius (mm) | fitted |
| `fr0dev` | relative SD of the maximum radius | 0.2 |
| `frimm` | factor on the maximum radius after t_imm | = `fvimm` |
| `frimmdev` | relative SD of the maximum radius after t_imm | 0.2 |
| `b` | inverse breathing amplitude (mm^-1): how far new lesions appear | fitted |
| `nata` | natality, rate of reinfection (day^-1 mm^-1) | fitted |
| `f` | factor on the wall distance that limits the radius | 1 |
| `tmax` | simulated days | 80 |
| `dt` | time step (days) | 0.1 |
| `n`, `a` | time course of the reinfection rate | 1.63, 0.035 |
| `rmin` | initial radius (mm) | 0.07 |
| `n0` | mean number of initial lesions | fitted |
| `fndev` | relative SD of the number of initial lesions | 0.2 |
| `tc` | day when the immune response starts | fitted |
| `ftcdev` | relative SD of t_imm | 0 |
| `tmin` | age at which lesions start growing (days) | 7 |

## Outputs (`results/`)

- `results.txt` - one line per day k = 0 ... tmax-1 (state at the end of day k + 1),
  comma-separated: k, then mean and SD over the simulations of: all lesions, NC, NM, VFC, NMD,
  NMI, mean main axis, MAC (columns 2-17), followed by further diagnostics not used in the paper
  (fusions, reinfections inside / outside the lobe, per-lobe counts; columns 18-51).
- `error.txt` - one line per run: E_NC, E_NM, E_NMI, E_NMD, E_VFC, E_A, E_B (errors against
  `config/experimental_results.txt` at days 21, 53, 77, the last weighted 5x;
  E_A = E_NC + E_NM + E_VFC, E_B = E_NC + E_NMI + E_NMD + E_VFC). If a simulation exceeds the
  model's limits (lesion counts, 600 s), the line is `100, ..., 100`.
- `prcc.txt` - one line per run: the 19 parameters, the daily means of NC, NM, NMD, NMI, VFC, MAC,
  and the errors.

Lung volume (319 cm^3) and trachea diameter (8.1 mm) of the virtual lung are constants in
`include/main.h`.
