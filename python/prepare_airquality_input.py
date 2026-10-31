"""Prepare an independently downloaded original Air Quality XLSX locally."""
import argparse
from copa2026.figures import prepare_airquality_input


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xlsx_path")
    parser.add_argument("output_dir", help="External local folder; existing inputs are never replaced")
    args = parser.parse_args(argv)
    path = prepare_airquality_input(args.xlsx_path, args.output_dir)
    print(path)
    return path


if __name__ == "__main__":
    main()
