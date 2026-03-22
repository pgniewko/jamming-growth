[![DOI](https://zenodo.org/badge/138428274.svg)](https://zenodo.org/badge/latestdoi/138428274)

# jamming-growth

2D numerical simulations of growing budding-yeast packings, derived shear-response
calculations, and `B_ext` probes.

The repository now uses a single growth code:

- `src/jamming_by_growth.f` generates the jammed packing at `phi_j`, the final
  packing at `phi_j + dphi`, and the lineage sidecar files for both endpoints
- `src/shear_yeast_linearshear.f` computes the shear-response data
- `src/box_compress_bext.f` computes the `B_ext` probe data

## Installation

Required tools:

- `bash`
- `make`
- `gfortran`
- Python 3.10+ with `venv`

If needed, install the system dependencies first.

macOS:

```bash
brew install gcc python
```

Ubuntu or Debian:

```bash
sudo apt-get update
sudo apt-get install -y gfortran python3 python3-venv make
```

From the repo root:

```bash
bash ./install.sh
```

That script:

- creates `.venv/`
- installs any Python packages listed in `requirements.txt` if there are any
- builds the Fortran binaries into `bin/`
- verifies Python imports only when third-party requirements are listed

If you want the installer to try installing missing system packages:

```bash
bash ./install.sh --install-system-deps
```

If you want to rebuild manually:

```bash
make
```

## Single-Job Examples

Run one growth job:

```bash
bash examples/run_jamming.sh --results-root output/example
```

By default that keeps the compressed endpoint packings and analysis files:

- `LF_JAMM_<base>.gz`
- `LF_DPHI_<base>.gz`
- `LINEAGE_LF_JAMM_<base>.gz`
- `LINEAGE_LF_DPHI_<base>.gz`
- `STATS_LF_JAMM_<base>`
- `STATS_LF_DPHI_<base>`
- `NC_<base>`
- `STEPLOG_<base>`
- `DIVLOG_<base>`
- `TRANSITIONS_<base>`
- `POSTJAMM_SUMMARY_<base>`
- `COHORT_INITIAL_FREE_<base>.gz`

`COHORT_INITIAL_FREE_<base>.gz` is a compact, one-row-per-bud summary for buds
that were unconstrained at jamming. It stores the observed `a*` event when a bud
becomes constrained, plus the final right-censored state for buds that remain
unconstrained at the saved `dphi`.

The raw growth trajectory is deleted by default. Keep it with:

```bash
bash examples/run_jamming.sh --results-root output/example --save-all-data
```

Run one shear job from an existing growth endpoint:

```bash
bash examples/run_shear.sh --results-root output/example
```

By default this keeps only `G_data_LF_DPHI_<base>` and the stdout log. Keep the
shear trajectory with:

```bash
bash examples/run_shear.sh --results-root output/example --save-all-data
```

## Sweep Launchers

The canonical sweep launcher is:

```bash
python3 scripts/run_pipeline.py
```

That pipeline runs, for each parameter tuple:

1. growth
2. shear
3. `B_ext`

Shared options:

```bash
python3 scripts/run_pipeline.py --n-cpus 8
python3 scripts/run_pipeline.py --force
python3 scripts/run_pipeline.py --save-all-data
```

The stage-specific sweep launchers are:

```bash
python3 scripts/run_growth.py
python3 scripts/run_shear.py
python3 scripts/run_bext.py
```

`run_shear.py` and `run_bext.py` operate only on tuples whose growth outputs
already validate successfully.

`run_bext.py` accepts a probe compression override:

```bash
python3 scripts/run_bext.py --dphi-probe 1e-8
```

## Monitoring

Monitor the current dataset with:

```bash
python3 scripts/monitor_completed_jobs.py
```

The monitor validates content, not just filenames. A tuple is reported as:

- growth-complete only if the compressed packings at `phi_j` and `phi_j + dphi`
  both exist and parse successfully, together with the required stats and
  lineage files
- shear-complete only if `G_data` exists and parses successfully
- `B_ext`-complete only if the `B_ext` data file exists and parses successfully
- pipeline-complete only if all three stages validate

## Restart Behavior

Without `--force`, the single-job wrappers and sweep launchers skip only outputs
that validate successfully.

If a tuple is incomplete, empty, truncated, or otherwise invalid, only that
tuple is cleaned and rerun:

- invalid growth: rerun growth and clear that tuple's derived shear and `B_ext`
  outputs first
- invalid shear only: rerun only shear for that tuple
- invalid `B_ext` only: rerun only `B_ext` for that tuple

With `--force`, the requested stage scope is rerun even if outputs already
validate.

## Output Layout

The standard output tree is:

- `output/growth/`
- `output/shear/`
- `output/bext/`
- `output/logs/growth/`
- `output/logs/shear/`
- `output/logs/bext/`

Completed packings are stored compressed. Optional large outputs retained with
`--save-all-data` are also stored compressed.

## Repo Layout

- `src/`: Fortran simulation sources
- `examples/`: single-job shell wrappers
- `scripts/`: sweep launchers and shared pipeline helpers
- `bin/`: compiled executables created by `make`
- `requirements.txt`: third-party Python requirements, currently empty

## References

1. [Biomechanical feedback strengthens jammed cellular packings](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.208102), P. Gniewek, C. Schreck, O. Hallatschek, Phys. Rev. Lett. (2019)
2. [Self-driven jamming in growing microbial populations](https://www.nature.com/articles/nphys3741), M. Delarue, J. Hartung, C. Schreck, P. Gniewek, L. Herminghaus, O. Hallatschek, Nature Phys. (2016)
