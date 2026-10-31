"""Run the prepared Air Quality stream and generate Figure 3."""

import os
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

import argparse
import json
from copa2026.figures import run_real


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", help="External local folder containing prepared airquality.npz")
    parser.add_argument("--output-root")
    args = parser.parse_args(argv)
    _, summary, directory = run_real(data_dir=args.data_dir, output_root=args.output_root)
    print(json.dumps(summary, indent=2))
    print(directory)
    return directory


if __name__ == "__main__":
    main()
