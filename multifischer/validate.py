"""Validation and normalisation tools."""

from pathlib import Path
import numpy as np
import re

from . import utils
from . import outlog
from . import columns as cols


WAVELENGTH_ALIASES = (
    'wavelength',
    'wavelength (nm)',
    'wavelength(nm)',
    'wavelength / nm',
    'wavelength/nm',
    'wavelength [nm]',
    'wavelength[nm]',
    'wl',
    'wl (nm)',
    'wl(nm)',
    'wl / nm',
    'wl/nm',
    'wl [nm]',
    'wl[nm]',
    'nm',
    '(nm)',
    '[nm]',
)

DARK_ALIASES = (
    'dark',
    'ambient',
    'stable',
    'e',
    'pure e',
    'z',
    'pure z',
    'pure',
    'nonirradiated',
    'non-irradiated',
    'not irradiated',
    'initial',
)


def args(parser, args):
    """Validate and normalise parsed command-line arguments in place.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser used to report command-line validation errors.
    args : argparse.Namespace
        Parsed arguments to validate and normalise.
    """

    if not args.input.is_file():
        parser.error(f'Input file not found: {args.input}')

    if args.out is None:
        args.out = Path(args.input.stem)
    else:
        args.out = Path(args.out)
    while args.out.is_dir():
        overwrite = input(f'Output folder "{args.out}" already exists. Press ENTER to overwrite, or enter a new output folder directory: ')
        if not overwrite:
            break
        args.out = Path(overwrite)

    args.qyratio_range = require_ordered_range(
        parser, '--qyratio_range', args.qyratio_range, positive=True
    )

    args.darkmax_range = require_ordered_range(
        parser, '--darkmax_range', args.darkmax_range, positive=True
    )

    args.meta_range = require_ordered_range(
        parser, '--meta_range', args.meta_range, positive=True
    )

    args.pss_range = require_ordered_range(
        parser, '--pss_range', args.pss_range, positive=False
    )

    args.uvvis_range = require_ordered_range(
        parser, '--uvvis_range', args.uvvis_range, positive=True
    )

    if not (0.0 < float(args.hdi) <= 1.0):
        parser.error('--hdi must be between 0 and 1.')

    if not np.isfinite(args.hdi_sigma) or float(args.hdi_sigma) <= 0.0:
        parser.error('--hdi-sigma must be a finite positive value.')
    
    if not np.isfinite(args.meta_min):
        parser.error('--meta-min must be finite.')

    args.metastate = str(args.metastate).strip()
    if not args.metastate:
        parser.error('--metastate must not be empty.')

    if args.singlepair:
        normalised_singlepair = []

        for value in args.singlepair:
            try:
                wavelength = utils.wavelength_value_to_float(value)
            except ValueError:
                parser.error('--singlepair values must be irradiation wavelengths, e.g. "365" or "365 nm".')
            if wavelength <= 0.0:
                parser.error('--singlepair wavelengths must be positive values.')

            normalised_singlepair.append(utils.wavelength_label(wavelength))

        args.singlepair = tuple(normalised_singlepair)
    
    if args.violin_irrs is not None:
        if any(str(value).strip().lower() == 'all' for value in args.violin_irrs):
            if len(args.violin_irrs) != 1:
                parser.error(
                    '--violin-irrs "all" cannot be combined with individual wavelengths.'
                )
            args.violin_irrs = 'all'
        
        else:
            normalised_violin_irrs = []

            for value in args.violin_irrs:
                try:
                    wavelength = utils.wavelength_value_to_float(value)
                except ValueError:
                    parser.error(
                        '--violin-irrs values must be irradiation wavelengths or "all".'
                    )
                
                label = utils.wavelength_label(wavelength)
                if label not in normalised_violin_irrs:
                    normalised_violin_irrs.append(label)
                
            args.violin_irrs = normalised_violin_irrs



def require_ordered_range(parser, arg_name, arg_values, *, positive=False):
    """Validate and return an ordered two-value numeric range.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser used to report command-line validation errors.
    arg_name : str
        Argument name used in validation messages.
    arg_values : iterable of float
        Proposed lower and upper limits.
    positive : bool, optional
        Require both limits to be greater than zero.

    Returns
    -------
    tuple[float, float]
        Validated lower and upper limits.
    """

    arg_lower, arg_upper = [float(arg_value) for arg_value in arg_values]

    if not np.all(np.isfinite([arg_lower, arg_upper])):
        parser.error(f'{arg_name} values must both be finite.')
    if positive and (arg_lower <= 0.0 or arg_upper <= 0.0):
        parser.error(f'{arg_name} values must both be greater than zero.')
    if arg_lower >= arg_upper:
        parser.error(f'{arg_name} must be supplied as (lower, upper), with lower < upper.')

    return arg_lower, arg_upper


def numeric_input(df, irr_wls):
    """Validate converted numeric spectral data.

    Parameters
    ----------
    df : pandas.DataFrame
        Normalised spectral data after numeric conversion.
    irr_wls : list[str]
        Irradiation-spectrum column names.

    Raises
    ------
    ValueError
        If wavelengths are missing, duplicated, non-finite, or insufficient, or
        if any absorbance spectrum contains missing, non-numeric, or infinite
        values.
    """

    if df[cols.WAVELENGTH].isna().any():
        raise ValueError(
            'The wavelength column contains missing or non-numeric values.'
        )
    if df[cols.WAVELENGTH].duplicated().any():
        raise ValueError('The wavelength column contains duplicate values.')
    if not np.all(np.isfinite(df[cols.WAVELENGTH].to_numpy(dtype=float))):
        raise ValueError('The wavelength column must contain only finite values.')
    if len(df) < 2:
        raise ValueError(
            'The input data must contain at least two wavelength rows.'
        )       
    if df[cols.DARK].isna().any():
        raise ValueError('The dark spectrum contains missing or non-numeric absorbance values.')
    if not np.all(np.isfinite(df[cols.DARK].to_numpy(dtype=float))):
        raise ValueError('The dark spectrum must contain only finite values.')

    empty_irrs = [col for col in irr_wls if df[col].isna().any()]
    if empty_irrs:
        raise ValueError(
            'The following irradiation spectra contain missing or non-numeric absorbance values: '
            f'{empty_irrs}.'
        )

    inf_irrs = [col for col in irr_wls if not np.all(np.isfinite(df[col].to_numpy(dtype=float)))]
    if inf_irrs:
        raise ValueError(
            f'The following irradiation spectra contain infinite values: {inf_irrs}.'
        )


def col_header(df, col_name, col_alias):
    """Detect and standardise one required input column.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data with stripped, lower-case column names.
    col_name : str
        Standard column name.
    col_alias : iterable of str
        Accepted aliases for the required column.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the detected column renamed to ``col_name``.

    Raises
    ------
    ValueError
        If no matching column or more than one matching column is found.
    """

    matching_colnames = [col for col in df.columns if col in col_alias]
    if len(matching_colnames) == 0:
        raise ValueError(f'Missing {col_name} data. Ensure input CSV contains a "{col_name}" column.')
    if len(matching_colnames) > 1:
        raise ValueError(f'Multiple columns matching "{col_name}" found. Remove duplicate columns.')

    found_col = matching_colnames[0]
    outlog.col_found(col_name, found_col)
    return df.rename(columns={found_col: col_name})


def irr_header(df):
    """Detect and normalise irradiation-spectrum column headings.

    Irradiation columns are identified from numeric headings with an optional
    ``nm`` suffix. Their headings are normalised to compact numeric labels and
    sorted numerically.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data containing standardised wavelength and dark-spectrum
        columns.

    Returns
    -------
    tuple[pandas.DataFrame, list[str]]
        DataFrame with normalised irradiation columns and numerically sorted
        irradiation wavelength labels.

    Raises
    ------
    ValueError
        If fewer than two irradiation columns are detected or normalisation
        produces duplicate wavelength labels.
    """

    irr_colnames = []
    rename_map = {}
    pattern = re.compile(r'^\s*[-+]?(?:\d+\.?\d*|\.\d+)\s*(?:nm)?\s*$', re.I)

    for col in df.columns:
        if col in (cols.WAVELENGTH, cols.DARK):
            continue
        elif pattern.match(str(col)):
            label = utils.wavelength_label(col)
            rename_map[col] = label
            irr_colnames.append(label)

    if len(irr_colnames) <= 1:
        raise ValueError(
            'Insufficient irradiation columns found. '
            f'Irradiation columns detected: {irr_colnames}. '
            'At least two irradiated PSS spectra are required in addition '
            'to the dark-state spectrum.'
        )

    if len(set(irr_colnames)) != len(irr_colnames):
        duplicates = sorted(
            {col for col in irr_colnames if irr_colnames.count(col) > 1}
        )
        raise ValueError(f'Duplicate irradiation columns found: {duplicates}.')

    df = df.rename(columns=rename_map)
    irr_colnames = sorted(irr_colnames, key=utils.wavelength_value_to_float)
    outlog.irrs_found(irr_colnames)
    return df, irr_colnames


def range_overlap_wl(wl_data, requested_range, range_name):
    """Require a requested wavelength range to contain experimental data."""

    lower, upper = requested_range

    overlaps = (
        (wl_data >= lower)
        & (wl_data <= upper)
    ).any()

    if not overlaps:
        raise ValueError(f'{range_name} ({lower}-{upper} nm) does not overlap with the experimental wavelength range.')


def singlepair(irr_wls, singlepair):
    """Require both requested single-pair wavelengths to exist in the data."""

    missing = [wl for wl in singlepair if wl not in irr_wls]

    if missing:
        raise ValueError(
            "The following --singlepair wavelength(s) were not found "
            f"in the input data: {', '.join(missing)}."
        )


def irradiation_wavelength_range(wl_data, irr_wls):
    """Require all irradiation wavelengths to lie within the measured spectrum."""

    spectral_min = wl_data.min()
    spectral_max = wl_data.max()

    outside_range = [
        irr_wl
        for irr_wl in irr_wls
        if not spectral_min <= float(irr_wl) <= spectral_max
    ]

    if outside_range:
        raise ValueError(
            'The following irradiation wavelengths lie outside the measured '
            f'spectral range ({spectral_min:g}-{spectral_max:g} nm): '
            f'{", ".join(outside_range)} nm.'
        )


def violin_irrs(irr_wls, requested_irrs):
    """Validate requested violin plot irradiation wavelengths exist."""

    if requested_irrs in (None, 'all'):
        return
    
    missing = [irr for irr in requested_irrs if irr not in irr_wls]

    if missing:
        raise ValueError(
            'The following --violin-irrs wavelengths were not found in the '
            f'input data: {missing}. Available wavelengths: {irr_wls}.'
        )