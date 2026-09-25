"""Run a list of CLASS models in parallel and save each one."""

import os
from concurrent.futures import ProcessPoolExecutor

from run_class import run, save


def _run_one(job):
    name, params, config, directory, with_rates = job
    try:
        out, meta = run(params, config, with_rates=with_rates)
    except Exception as e:  # classy raises CosmoSevereError / CosmoComputationError
        return name, str(e).strip().splitlines()[-1]
    save(name, out, meta, directory)
    return name, None


def run_grid(jobs, config, directory, with_rates=True, workers=None):
    """jobs: list of (name, params). Returns {name: error message or None}.

    CLASS is itself OpenMP-parallel, so the process count is kept small by default.
    """
    workers = workers or max(1, (os.cpu_count() or 4) // 4)
    tasks = [(name, params, config, directory, with_rates) for name, params in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return dict(pool.map(_run_one, tasks))
