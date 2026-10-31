"""Run or resume the joint Figure 1/2 sensitivity experiment."""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
from copa2026.figures import run_full_figure_sweeps


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir")
    parser.add_argument("--datasets", nargs="+", default=["protein", "particle", "sensorless", "credit"])
    parser.add_argument("--k-list", type=int, nargs="+", default=list(range(1, 101)))
    parser.add_argument("--lambda-list", type=float, nargs="+", default=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    parser.add_argument("--no-distance", action="store_true")
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--epsilon", type=float, default=100.)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--output-root", help="New run directory or existing matching run to resume")
    parser.add_argument("--max-jobs", type=int, default=None, help="Limit new jobs in this invocation")
    args = parser.parse_args(argv)
    rows, directory = run_full_figure_sweeps(data_dir=args.data_dir, dataset_names=args.datasets,
        k_list=args.k_list, lambda_list=args.lambda_list, include_distance=not args.no_distance,
        reps=args.reps, epsilon=args.epsilon, workers=args.workers, output_root=args.output_root,
        max_jobs=float("inf") if args.max_jobs is None else args.max_jobs)
    print(f"Completed checkpoints: {len(rows)}\n{directory}")
    return directory


if __name__ == "__main__":
    main()
