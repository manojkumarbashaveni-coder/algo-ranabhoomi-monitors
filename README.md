# Monitor Placement — ALGO RANABHOOMI submission

This repository contains my submission for the **MONITORS** challenge.  The
only solution implementation is `adapters/mine.py`, which defines the required
`MySolver` class.  The surrounding public harness files are unmodified copies
of the challenge handout so the solution can be inspected and run locally.

## Problem

Choose exactly `k` permitted monitor sites.  Each demand point has a population
weight and pays its weighted distance to the nearest selected site.  Distances
outside the coverage radius are charged a fixed penalty instead.  Candidate
sites are not necessarily demand points.

The benchmark preserves the best valid plan at 5%, 20%, 50%, and 100% of a
five-second budget, so both early and final solution quality matter.

## Algorithm

The solver is an anytime variable-neighborhood local search:

1. It submits a deterministic valid fallback plan immediately.
2. It greedily adds sites, selecting the candidate with the best reduction in
   weighted charged distance.
3. It repeatedly applies the best exact one-site exchange until reaching a
   one-swap local optimum.
4. It removes two chosen sites and greedily repairs with two different
   replacements.  All pairs are considered for small `k`; larger instances use
   a deterministic sampled set to respect the time budget.
5. It uses bounded 2–4-site destroy-and-greedy-repair restarts to escape local
   optima, always retaining and submitting the best valid plan found.

The two-site neighborhood helps when two monitor choices must change together:
a single swap can be temporarily worse even though the paired replacement is
better overall.

## Results

| Version | Public score |
| --- | ---: |
| Published greedy + one swap baseline | 55/100 |
| Greedy + repeated swaps + random restart | ~86/100 |
| Final solver | ~93/100 |

The final public-equivalent check measured 92.997/100.  Results on the
organizer's Linux evaluation environment may vary slightly with timing.

## Run

Use Python 3.  The supplied runner uses Unix process-group APIs, so run it on
Linux, macOS, or WSL rather than native Windows.

```bash
python self_check.py --adapter adapters.mine:MySolver
python run.py --adapter adapters.mine:MySolver --out report.json
```

`report.json` is generated locally and is intentionally ignored by Git.

## Repository contents

- `adapters/mine.py` — final submitted solver.
- `adapter.py`, `data.py`, `validator.py` — unchanged challenge interfaces and
  instance/validation logic.
- `run.py`, `self_check.py`, `metrics.py`, `_paths.py`, `public_reference.json`,
  and `benchkit/` — unchanged public harness required to run the commands above.
