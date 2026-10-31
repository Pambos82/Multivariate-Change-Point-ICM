"""Render complete Figures 1 and 2 from a finished Performance.csv grid."""
import argparse
from copa2026.figures import reconstruct_copa_figures


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("performance_input", help="Performance.csv or its containing directory")
    parser.add_argument("--export-pdf", action="store_true")
    args = parser.parse_args(argv)
    directory = reconstruct_copa_figures(args.performance_input, export_pdf=args.export_pdf)
    print(directory)
    return directory


if __name__ == "__main__":
    main()
