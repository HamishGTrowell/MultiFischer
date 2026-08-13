"""Command-line interface for xFischer PSS analysis."""

import sys
import argparse
from pathlib import Path

from . import workflow
from . import validate
from . import outlog


def parse_args(argv):
    """Parse and validate command-line arguments.

    Parameters
    ----------
    argv : iterable of str
        Command-line arguments excluding the program name.

    Returns
    -------
    argparse.Namespace
        Validated and normalised command-line arguments.

    Raises
    ------
    SystemExit
        If parsing or validation fails, or if help is requested.
    """

    parser = argparse.ArgumentParser(
        description='xFischer: Extended Fischer Method PSS Analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.suggest_on_error = True

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        metavar='INPUT.CSV',
        help=(
            'Required: Input CSV file directory for UV-Vis data with columns: "Wavelength", "Dark", and '
            'irradiation spectra labelled by wavelength, e.g. "365" or "365 nm".'
        ),
    )

    parser.add_argument(
        '--out',
        type=Path,
        metavar='OUT/',
        help='Optional: Output folder directory for saving results. Matches input filename if not specified.'
    )

    parser.add_argument(
        '--singlepair',
        dest='singlepair',
        nargs=2,
        default=False,
        metavar=('IRR1', 'IRR2'),
        help=(
            'Optional: Irradiation wavelengths for analysis of a single irradiation pair, '
            'e.g. --singlepair 340 365.'
        ),
    )

    parser.add_argument(
        '--qyratio-range',
        dest='qyratio_range',
        type=float,
        nargs=2,
        default=(0.5, 2.0),
        metavar=('MIN_X', 'MAX_X'),
        help='Optional: Quantum-yield-ratio range (X) for sensitivity analysis.',
    )

    parser.add_argument(
        '--darkmax-range',
        dest='darkmax_range',
        type=float,
        nargs=2,
        default=(275.0, 800.0),
        metavar=('MIN_NM', 'MAX_NM'),
        help='Optional: Wavelength range (nm) used to locate lambda_max(Dark), selected by maximum absorbance within specified range.',
    )

    parser.add_argument(
        '--meta-range',
        dest='meta_range',
        type=float,
        nargs=2,
        default=(275.0, 600.0),
        metavar=('MIN_NM', 'MAX_NM'),
        help=(
            'Optional: Wavelength range (nm) over which the extrapolated metastable-state absorbance '
            'must be greater than meta_min. Used to filter irradiation pairs which lead to unphysical negative absorbance of the extrapolated metastate.'
        ),
    )

    parser.add_argument(
        '--meta-min',
        dest='meta_min',
        type=float,
        default=-0.05,
        metavar='MIN_ABS',
        help=(
            'Optional: Minimum allowed extrapolated metastable-state absorbance '
            'in meta_range. Used to filter irradiation pairs which lead to unphysical negative absorbance of the extrapolated metastate.'
        ),
    )

    parser.add_argument(
        '--pss-range',
        dest='pss_range',
        type=float,
        nargs=2,
        default=(-0.05, 1.05),
        metavar=('MIN_PSS', 'MAX_PSS'),
        help='Optional: Accepted PSS range (as a decimal value between 0 and 1). Used to filter irradiation pairs which lead to unphysical values for the PSS far outside 0-1.',
    )

    parser.add_argument(
        '--hdi',
        type=float,
        default=0.9,
        metavar='HDI',
        help='Optional: Fraction of candidate PSS datapoints used to define the highest density interval (HDI) core. Used to filter irradiation pairs with anomalous results.',
    )

    parser.add_argument(
        '--hdi-sigma',
        dest='hdi_sigma',
        type=float,
        default=5.0,
        metavar='HDI_SIGMA',
        help='Optional: Multiple of HDI-core standard deviations used to define extended HDI bounds. Used to filter irradiation pairs with anomalous results.',
    )

    parser.add_argument(
        '--uvvis-range',
        dest='uvvis_range',
        type=float,
        nargs=2,
        default=(200.0, 600.0),
        metavar=('MIN_NM', 'MAX_NM'),
        help='Optional: Wavelength range (nm) for the UV-Vis plot.',
    )

    parser.add_argument(
        '--violin-irrs',
        dest='violin_irrs',
        nargs='+',
        default=None,
        metavar='IRR',
        help=(
            'Optional: Irradiation wavelengths shown in the PSS distribution '
            'violin plot, e.g. --violin-irrs 340 365, or "all". If omitted, '
            'the irradiation spectrum with the highest PSS is used.'
        ),
    )

    parser.add_argument(
        '--sensitivity',
        action='store_true',
        help='Optional: Include sensitivity ranges in final UV-Vis and PSS plots.',
    )

    parser.add_argument(
        '--metastate',
        default='Z',
        metavar='METASTATE',
        help='Optional: Identity of the metastable isomer, e.g. "Z" for azobenzene.',
    )

    args = parser.parse_args(list(argv))

    validate.args(parser, args)

    return args


def main():
    """Command-line entry point.
    Run the command-line workflow and exit with its return code.

    Raises
    ------
    SystemExit
        Raised with 0 after a successful run or 1 after a handled analysis error.    
    """

    args = parse_args(sys.argv[1:])
    args.out.mkdir(parents=True, exist_ok=True)
    outlog.setup(args.out)

    try:
        return_code = workflow.run(args)

    except ValueError as exc:
        outlog.errorln(str(exc))
        outlog.error_termination()
        print(f'ERROR: {exc}', file=sys.stderr)
        return_code = 1

    finally:
        outlog.close()

    raise SystemExit(return_code)


if __name__ == '__main__':
    main()
