"""Run the four-dataset benchmark or replay saved benchmark p-values."""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse

from copa2026.benchmarks import run_benchmarks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--epsilon", type=float, default=100)
    parser.add_argument("--threshold", type=float, default=100)
    parser.add_argument("--data-dir")
    parser.add_argument("--output-root")
    parser.add_argument("--replay-file")
    args = parser.parse_args(argv)
    return run_benchmarks(reps=args.reps, epsilon=args.epsilon, threshold=args.threshold,
                          data_dir=args.data_dir, output_root=args.output_root,
                          replay_file=args.replay_file)


if __name__ == "__main__":
    main()
