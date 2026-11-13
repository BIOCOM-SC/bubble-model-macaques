# Bubble Model - TB lesions in rhesus macaques

Data and code of the paper

> F. Muñoz, M. Català, I. Nogueira, B. Eguzkitza, S. Sharpe, C. Prats, P.-J. Cardona.
> *Assessment of protective mechanisms in the evolution of* Mycobacterium tuberculosis *infection
> in rhesus macaques using a parameterized computational model.* Scientific Reports (2026).
> DOI: _to be added_

The Bubble Model (BM) simulates the growth, merging and spread of tuberculosis lesions in a virtual
macaque lung. Its parameters are fitted to CT data of rhesus macaques infected with
*M. tuberculosis* (Nogueira et al. 2023), and the fitted model is then used for in-silico
experiments, sensitivity analyses and a comparison with BCG-vaccinated animals.

## Contents

- [`data/`](data/README.md): experimental data (CT lesion measurements, lung volumes) and the
  script that turns them into the model's targets (Table 2).
- [`lung/`](lung/README.md): the virtual macaque lung (lobe and segment surfaces, terminals of the
  bronchial tree and the distances between them).
- [`model/`](model/README.md): the Bubble Model (C, MPI + OpenMP).
- [`runs/`](runs/README.md): configuration and launcher of every simulation run of the paper (fits,
  experiments, sweeps, PRCC), for a Slurm cluster.
- [`analysis/`](analysis/README.md): collection of the run outputs into tables, and the gnuplot
  scripts of every figure.

The pipeline:

```
data/raw/ ── process.py ──> data/processed/ (targets) ─┐
lung/macaque/ (virtual lung) ──────────────────────────┼──> runs/launch.py ──> runs/work/
model/ (C) ────────────────────────────────────────────┘                            │
                                                                           analysis/collect.py
                                                                                    │
analysis/out/ (figures) <── plots/*.plt ── analysis/data/, analysis/tables/ <───────┘
```

The run outputs themselves (`runs/work/`, ~1 GB) are not included. Everything the figures and
tables need is: `analysis/data/` and `analysis/tables/` hold the collected results of the runs.

## How to use it

**1. Redraw the figures** from the included results (gnuplot only):

```sh
cd analysis
./plot_all.sh --plots            # all figures -> analysis/out/ (EPS, PDF, PNG)
./plot_all.sh --plots fits       # one figure
```

**2. Recompute the experimental targets** from the raw data (Python):

```sh
python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
python data/process.py                 # unvaccinated macaques -> data/processed/
python data/process.py --cohort RIDV   # intradermal-BCG-vaccinated macaques -> data/processed/ridv/
```

**3. Rerun the simulations.** Unpack the lung distance file once, build the model, and run:

```sh
xz -dk lung/macaque/D.bin.xz           # 60 MB -> 271 MB
make -C model                          # local build (see model/README.md)
python runs/launch.py run experiments --test --local   # small test on a laptop (seconds)
```

The runs of the paper need a cluster: about 3.5 x 10^8 simulated animals, about 40 hours on two
112-core nodes in total, split into four jobs that can run in parallel (the longest takes ~15 h).
On a Slurm cluster, copy `hpc_site.example.toml` to `hpc_site.toml`, fill it in, and submit
everything with `python runs/launch.py submit all` (see [runs/README.md](runs/README.md)).
Then `cd analysis && ./plot_all.sh` collects the new outputs and redraws all figures.

## Requirements

- Python >= 3.11 with the packages in `requirements.txt` (data processing and result collection;
  the run launcher itself uses the standard library only)
- gnuplot >= 5.4, `epstopdf` and Ghostscript (figures)
- an MPI C compiler with OpenMP and C23 support (model), `xz`
- [Dakota](https://dakota.sandia.gov) >= 6 (parameter fits only)

## Data sources

- CT lesion data: Nogueira, I. *et al.* Pathogens **12**, 236 (2023),
  [doi:10.3390/pathogens12020236](https://doi.org/10.3390/pathogens12020236), Supplementary File S2,
  redistributed unchanged under its CC BY 4.0 licence.
- Lung volumes (one per animal, CT at week 3): provided by the authors of Nogueira et al. (2023).
- Virtual lung: described in the Methods of the paper.

## Licence

Code: MIT (see `LICENSE`). Third-party data in `data/raw/`: CC BY 4.0, as stated above.

## Citation

Please cite the paper above. 
