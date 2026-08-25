# Is the lobe order robust to the coarse-graining?

Internal robustness check, decoupled from the report. The chapter derives the
uniform-rate activation order of the four lobes, temporal, occipital,
parietal, frontal, on the 83-region graph built by filtering the fine
Budapest graph at five occurrences and aggregating the 1015 parcels by
FreeSurfer parent region. This study repeats the reference tau run on the
fine graph itself, with no aggregation and nothing retuned, and aggregates
the lobes only after the simulation.

## Setup

`study-fine-graph-order.py`, run from the repository root. Same model as the
reference run of benchmark 27: metric-graph FEM with the lumped mass, one P1
element per unit-length connection, backward Euler with Newton to 1e-11,
alpha = 0.5 uniform, rho = 0.005, dt = 0.4 over 80 years, seed c0 = 0.1 on
the parcels whose parent region is an entorhinal cortex (5 parcels), lobes
assigned to parcels through the parent-region names with the classify rule
of the report. The integrator is validated first: on the 83-region graph it
reproduces the four stored lobe crossings to 4e-15 years.

## Result

```text
lobe crossings of the mean curves [years]
  lobe       fine, parcel mean   fine, region mean   coarse (stored)
  temporal        29.04               27.15               16.50
  occipital       40.70               41.77               21.65
  parietal        38.77               39.01               22.95
  frontal         41.19               41.37               25.87

orders
  coarse                 temporal, occipital, parietal, frontal   (span  9.37)
  fine, parcel mean      temporal, parietal, occipital, frontal   (span 12.15)
  fine, region mean      temporal, parietal, frontal, occipital   (span 14.62,
                         frontal and occipital 0.40 years apart, near a tie)

per-region activation times, fine against coarse: Spearman 0.635 over 83
```

Robust across the two levels: the temporal lobe is first, the frontal lobe
is never early (last, or third within 0.4 years of last) and the lobes
separate by years at this scaling. Sensitive to the coarse-graining: the
position of the occipital lobe, second on the aggregated graph, third or
last on the fine one, so the specific middle of the sequence
temporal-occipital-parietal-frontal is a property of the aggregated
representation and not of the fine connectivity.

## Where the difference is born

- Not in the lobe-to-lobe totals: aggregation sums the same fine weights, so
  the six totals are identical at the two levels to the printed digit
  (74.99, 55.58, 1.03, 190.46, 98.40, 2.27).
- Not in the direct couplings either: per parcel, the occipital lobe is
  still the better-connected one to the temporal lobe (coupling over lumped
  mass 0.0130 mean against 0.0050 for the parietal), so the coarse
  explanation, occipital ahead because its share per unit is larger, points
  the same way at the fine level, yet the order flips there.
- The flip therefore comes from structure the aggregation removes or
  concentrates: the 7895 intra-region fine connections, which exist on the
  fine graph and are deleted by the aggregation, and the multi-hop inflow,
  which favours the parietal lobe (407 units of non-temporal inflow against
  108 for the occipital lobe, mostly through the frontal-parietal 190 and
  occipital-parietal 98 bundles) once the wave is no longer forced through
  single region-vertices.
- The seed also shrinks, 0.00025 of the metric mass against 0.00102 on the
  coarse graph, which delays everything (temporal 27-29 against 16.5 years)
  but does not by itself reorder lobes.

## Reading

The report's claim is stated for the reconstructed 83-region graph and this
check does not contradict it; what it bounds is its scope. Temporal first
and frontal late survive the change of representation; the exact position
of the occipital lobe does not, so the middle of the connectivity-driven
sequence should not be read as a property of the brain's fine connectivity.
No report text or figure uses this study.

## Files

- `study-fine-graph-order.py` — the study, self-validating.
- `results/fine_lobe_crossings.csv` — the four lobes at both aggregations.
- `results/fine_region_times.csv` — the 83 per-region times, fine and coarse.
