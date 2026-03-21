[![DOI](https://zenodo.org/badge/138428274.svg)](https://zenodo.org/badge/latestdoi/138428274)

# jamming-growth

2D numerical simulations of growing budding-yeast packings. The repository
contains two Fortran programs:

- `src/jamming_by_growth.f` to generate jammed packings during growth
- `src/shear_yeast_linearshear.f` to shear a saved packing

## Installation

If the required tools are not installed, use one of these commands first.

macOS:

```bash
brew install gcc python
```

Ubuntu or Debian:

```bash
sudo apt-get update
sudo apt-get install -y gfortran python3 python3-venv make
```

Required tools:

- `bash`
- `make`
- `gfortran`
- Python 3.10+ with `venv`

From the repo root:

```bash
bash ./install.sh
```

That script:

- creates `.venv/`
- installs Python packages from `requirements.txt`
- builds both Fortran binaries into `bin/`
- verifies that the required Python packages can be imported
- does not activate the environment

If you want the installer to try installing missing system packages:

```bash
bash ./install.sh --install-system-deps
```

If you want to rebuild the binaries manually after setup:

```bash
make
```

If you want to activate the environment after installation:

```bash
source .venv/bin/activate
```

## Running Examples

The repository includes two example wrappers in `examples/`.

Run a single growth simulation:

```bash
bash examples/run_jamming.sh --results-root output/example-growth
```

This creates the output and log directories automatically under the chosen
results root.
The growth run also writes growth-state files and diagnostics there, including
`LF_JAMM_*`, `LF_DPHI_*`, `NC_*`, and `STEPLOG_*`.

Run a shear simulation from the saved growth state produced by the growth run:

```bash
bash examples/run_shear.sh --results-root output/example-growth
```

This also creates its output and log directories automatically under the chosen
results root.

## Full growth + shear sweep

Run the full growth + shear sweep:

```bash
python3 scripts/run_growth_shear.py
```

Choose the number of concurrent jobs:

```bash
python3 scripts/run_growth_shear.py --n-cpus 8
```

Force reruns even when outputs already exist:

```bash
python3 scripts/run_growth_shear.py --force
```

Monitor completed growth + shear jobs:

```bash
python3 scripts/monitor_completed_jobs.py
```

Run the compact lineage-growth stage after building:

```bash
make
python3 scripts/run_lineage_growth.py --sizes 8,15 --p0s 1e-3,1e-4 --dphis 1e-3,1e-2 --seeds 1201,1202
```

Choose concurrency or force reruns:

```bash
python3 scripts/run_lineage_growth.py --sizes 8 --p0s 1e-3 --dphis 1e-3 --seeds 1201 --n-cpus 4 --force
```

Use the built-in large-system presets for the lineage pilot and the follow-on
paper sweep:

```bash
python3 scripts/run_lineage_growth.py --preset study-small10 --n-cpus 4
python3 scripts/run_lineage_growth.py --preset pilot-large --n-cpus 4
python3 scripts/run_lineage_growth.py --preset full-large --n-cpus 6
```

The lineage run writes the usual growth endpoints plus lineage-specific raw
data for post-jamming depletion analysis:

- `LINEAGE_LF_JAMM_*` and `LINEAGE_LF_DPHI_*` for endpoint lineage snapshots
- `DIVLOG_*` for division events
- `TRANSITIONS_*` for tracked initial-free-bud transitions
- `POSTJAMM_SUMMARY_*` for accepted post-jamming summary rows (`phi`, counts,
  `chi_c`, and tracked-reservoir depletion state)

Generate the lineage pilot review tables and figures from the raw lineage
outputs:

```bash
.venv/bin/python scripts/plot_lineage_support.py --preset study-small10
.venv/bin/python scripts/plot_lineage_support.py --preset pilot-large
```

Keep shear trajectory files after successful jobs:

```bash
python3 scripts/run_growth_shear.py --keep-all-output
```

## Standalone B_ext pass

Run the B_ext pass over existing growth packings:

```bash
python3 scripts/run_bext.py
```

Use the help output to inspect all available options, including concurrency,
reruns, probe size, and output retention:

```bash
python3 scripts/run_bext.py --help
```

Generate the manuscript-support figure set from the existing growth, shear,
and `B_ext` outputs:

```bash
.venv/bin/python scripts/plot_paper_support.py
```

By default this writes the current canonical figures and CSV summaries under
`output/plots/`.

## Repo Layout

- `src/`: Fortran simulation codes
- `examples/`: shell wrappers for single growth and shear runs
- `bin/`: compiled executables created by `make`
- `requirements.txt`: Python dependencies for current and future scripts

## References

1. [Biomechanical feedback strengthens jammed cellular packings](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.208102), P. Gniewek, C. Schreck, O. Hallatschek, Phys. Rev. Lett. (2019)
2. [Self-driven jamming in growing microbial populations](https://www.nature.com/articles/nphys3741), M. Delarue, J. Hartung, C. Schreck, P. Gniewek, L. Herminghaus, O. Hallatschek, Nature Phys. (2016)
