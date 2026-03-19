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

## Numerical experiments

Build the B_ext executable:

```bash
make box_compress_bext
```

Run the B_ext pass over existing growth packings:

```bash
python3 scripts/run_bext.py
```

Choose the number of concurrent B_ext jobs:

```bash
python3 scripts/run_bext.py --n-cpus 8
```

Force B_ext reruns even when valid raw outputs already exist:

```bash
python3 scripts/run_bext.py --force
```

Override the fixed B_ext probe compression:

```bash
python3 scripts/run_bext.py --dphi-probe 1e-5
```

## Repo Layout

- `src/`: Fortran simulation codes
- `examples/`: shell wrappers for single growth and shear runs
- `bin/`: compiled executables created by `make`
- `requirements.txt`: Python dependencies for current and future scripts

## References

1. [Biomechanical feedback strengthens jammed cellular packings](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.208102), P. Gniewek, C. Schreck, O. Hallatschek, Phys. Rev. Lett. (2019)
2. [Self-driven jamming in growing microbial populations](https://www.nature.com/articles/nphys3741), M. Delarue, J. Hartung, C. Schreck, P. Gniewek, L. Herminghaus, O. Hallatschek, Nature Phys. (2016)
