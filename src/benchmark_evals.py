"""Benchmark likelihood-evaluation throughput for different worker/thread splits.

Usage: uv run python src/benchmark_evals.py <likelihoods> <workers>x<threads> [...]
e.g.   uv run python src/benchmark_evals.py plik 45x4 90x2 30x6
Each configuration evaluates 3 x workers points near the Planck best fit and
reports evaluations per second (after one warm-up evaluation per worker).
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from likelihood_setup import info  # noqa: E402

_MODEL = None


def _init(threads, likelihoods):
    os.environ["OMP_NUM_THREADS"] = str(threads)
    global _MODEL
    from cobaya.model import get_model
    _MODEL = get_model(info(likelihoods, u=1e-4))


def _eval(seed):
    rng = np.random.default_rng(seed)
    p = {n: (v["ref"]["loc"] if isinstance(v.get("ref"), dict) else v.get("ref"))
         for n, v in _MODEL.parameterization.sampled_params_info().items()}
    p["omega_dm"] *= 1 + 1e-3 * rng.standard_normal()
    t = time.time()
    _MODEL.logposterior(p)
    return time.time() - t


def main():
    likelihoods = sys.argv[1]
    for spec in sys.argv[2:]:
        workers, threads = map(int, spec.split("x"))
        with ProcessPoolExecutor(workers, initializer=_init, initargs=(threads, likelihoods)) as pool:
            list(pool.map(_eval, range(workers)))               # warm-up
            t0 = time.time()
            times = list(pool.map(_eval, range(1000, 1000 + 3 * workers)))
            wall = time.time() - t0
        print(f"{spec:>6}: {3 * workers / wall:6.1f} evals/s  (single eval {np.median(times):.2f} s median)", flush=True)


if __name__ == "__main__":
    main()
