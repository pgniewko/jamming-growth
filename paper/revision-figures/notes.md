# Codex task: add supplementary post-jamming division-count figure

Goal: quantify **how many post-jamming cell divisions occur from \phi_J onward**, and how that activity is split between the interval **\phi_J \to \phi_2** and the interval **beyond \phi_2**.

Use only the already-saved lineage outputs. Do **not** add new simulations or modify the Fortran.

## Main observable
For each run, construct the cumulative count of **all** post-jamming divisions
\[
N_{\mathrm{div}}^{\mathrm{all}}(\Delta\phi)
=
\#\{\text{DIVLOG events with } \phi_{\mathrm{event}}-\phi_J \le \Delta\phi\},
\qquad
\Delta\phi = \phi-\phi_J .
\]

This is the **main curve**.

Plot it from `Delta phi = 0` (the jamming point) all the way to the last available post-jamming point.
Do **not** stop at `phi_2`.

## Secondary observable
Also compute the binned all-division count rate
\[
r_{\mathrm{div}}^{\mathrm{all}}(\Delta\phi)
=
\frac{\Delta N_{\mathrm{div}}^{\mathrm{all}}}{\Delta\phi},
\]
using fixed `Delta phi` bins.

This is a secondary panel only; it tells us **where** the division activity is concentrated.

## Reinjection subset
In addition, build the cumulative count for the direct reinjection subset:
\[
N_{\mathrm{div}}^{\mathrm{track}}(\Delta\phi)
=
\#\{\text{TRANSITIONS events with event\_code == 3 and } \phi-\phi_J \le \Delta\phi\}.
\]

This should be a smaller companion panel or inset.
It is the direct check on the initial-free-bud reinjection correction.

## Statistics to show
For fixed `(L, P0)`:

- For the cumulative counts, summarize over runs by **median** and **IQR**.
- For the binned rate, use the **pooled mean count per bin** with a **bootstrap confidence band over runs**.

The main y-axis should be the **count per run**.
If multiple `L` values are shown on the same axes, also save a normalized version
\[
N_{\mathrm{div}}^{\mathrm{all}} / N_J
\]
for cross-size comparison, but keep the raw count as the main supplement figure.

## Must report explicitly
For each run, compute:

- `N_div_all_at_phi2`
- `N_div_all_final`
- `N_div_all_after_phi2 = N_div_all_final - N_div_all_at_phi2`
- `frac_before_phi2 = N_div_all_at_phi2 / N_div_all_final` when `N_div_all_final > 0`

This is the key summary, because it tells us directly how much division activity happens:

1. between `phi_J` and `phi_2`
2. after `phi_2`

If a run does **not** exhaust the tracked reservoir by the endpoint, keep the cumulative curve to the endpoint, mark it as `exhausted_by_endpoint = 0`, and exclude it from the pre/post-`phi_2` split summary.

## Inputs already available
Use:

- `DIVLOG_*` for **all** post-jamming division events. :contentReference[oaicite:4]{index=4}
- `POSTJAMM_SUMMARY_*` for the post-jamming trajectory in `phi` and the tracked initial-free counts. 
- existing depletion logic to get the observed `phi_2` from the first step where `n_initial_free_active == 0`. :contentReference[oaicite:6]{index=6}
- `TRANSITIONS_*` with `event_code == 3` for the tracked-initial-free reinjection subset. :contentReference[oaicite:7]{index=7}

## Required code change
The current `load_divlog_summary()` only stores the **total number of DIVLOG rows per file**.
Replace or supplement it with an **event-level DIVLOG loader** that stores at least:

- run metadata
- `step`
- `phi`

Then compute
\[
\Delta\phi_{\mathrm{event}} = \phi_{\mathrm{event}} - \phi_J .
\]

`phi_J` should come from the first row of the matching `POSTJAMM_SUMMARY_*` run. :contentReference[oaicite:8]{index=8}

## Figure layout
Create `fig_lineage_division_counts` with three panels.

### Panel A — main panel
Median cumulative all-division count
\[
N_{\mathrm{div}}^{\mathrm{all}}(\Delta\phi)
\]
with IQR bands.

- x-axis: `Delta phi`
- y-axis: cumulative division count
- color by `P0`
- one subplot per `L`
- mark the median observed `Delta phi_2` for each `P0` with a color-matched vertical dashed line

### Panel B — where activity is concentrated
Binned all-division count rate
\[
r_{\mathrm{div}}^{\mathrm{all}}(\Delta\phi).
\]

- same x-axis range
- same `L` / `P0` organization
- extend beyond `phi_2` to the full available endpoint

### Panel C — before vs after `phi_2`
Bar plot (or point plot) of the median counts:

- divisions accumulated by `phi_2`
- divisions accumulated after `phi_2`

for each `P0` and `L`.

Optional inset:
\[
N_{\mathrm{div}}^{\mathrm{track}}(\Delta\phi)
\]
for the tracked initial-free subset only.

## Outputs
Write:

- `division_counts_cumulative.csv`
- `division_counts_rate.csv`
- `division_counts_summary.csv`

In the summary file include, for each `(L, P0, seed)`:

- `phi_j`
- `phi_2_obs` / `delta_phi_2_obs` if present
- `exhausted_by_endpoint`
- `N_div_all_at_phi2`
- `N_div_all_final`
- `N_div_all_after_phi2`
- `frac_before_phi2`
- `N_div_track_at_phi2`
- `N_div_track_final`

## Caption intent
State clearly:

This figure is a supplementary diagnostic of **when** post-jamming division activity becomes appreciable.
It does **not** by itself prove force-network breaking, but it identifies the `Delta phi` window where division-induced restructuring is most likely to matter mechanically, and it separately shows whether the direct reinjection subset remains small.
