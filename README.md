# class_baryons

CMB constraints on compact, weakly photon-coupled matter (e.g. dense baryonic clumps) using the [CLASS](https://github.com/lesgourg/class_public) Boltzmann code.

The central question: can a pressureless population of optically thick clumps stand in for some or all of the dark matter, and what minimum surface density must the clumps have for the CMB to remain CDM-like? The project also asks whether such macroscopic objects share the Thomson-like photon collision operator assumed in standard dark-matter–photon analyses.

## Documents

- [Experiment plan](doc/cmb_compact_baryon_clumps_class_experiment_plan.md)
- [Results summary](doc/RESULTS.md): current results and milestone status
- [Experimental record](doc/experiment_log/): dated, append-only run records
- [Implementation notes](doc/IMPLEMENTATION_NOTES.md): lessons for future implementers

## Quick start

```bash
scripts/setup_class.sh                 # clone CLASS v3.4.0 (pinned), build classy, uv sync
uv run python src/run_baseline.py      # M1: baseline LCDM
uv run python src/validate_m2a.py      # M2a: zero-coupling idm reproduces LCDM
uv run python src/validate_m2b.py      # M2b: cold-mass limit in m_idm
```
