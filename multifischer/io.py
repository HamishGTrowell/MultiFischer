"""Input/output helpers for loading UV-Vis data."""

import pandas as pd

from . import outlog
from . import validate
from . import columns as cols


def load_data(input_path, uvvis_range, darkmax_range, meta_range, singlepair, violin_irrs):
    """Load, validate, and normalise UV–Vis data from CSV.

    Required wavelength and dark-spectrum columns are detected from accepted
    aliases. Irradiation-spectrum columns are detected from wavelength-like
    headings and normalised to numeric string labels. All retained data are
    converted to numeric values and sorted by decreasing wavelength.

    Parameters
    ----------
    input_path : str or pathlib.Path
        Path to the input CSV file.
    uvvis_range : tuple[float, float]
        Wavelength range, in nm, used for UV-Vis plotting and data reporting.
    darkmax_range : tuple[float, float]
        Wavelength range, in nm, used to locate the dark-spectrum absorbance
        maximum.
    meta_range : tuple[float, float]
        Wavelength range, in nm, used for the metastate filter.
    singlepair : tuple[str, str] or False
        Requested single-pair irradiation wavelengths, or ``False`` for
        multipair analysis.
    violin_irrs : str or tuple[str] or ``all``
        Irradiation wavelengths to display violin plots for.

    Returns
    -------
    tuple[pandas.DataFrame, list[str]]
        Normalised spectral DataFrame and numerically sorted irradiation
        wavelength labels.

    Raises
    ------
    ValueError
        If the CSV cannot be loaded or its contents fail validation.
    """
    
    try:
        data_df = pd.read_csv(input_path)
    except Exception as exc:
        raise ValueError('Failed to load input data.') from exc
    outlog.logln('Input data loaded successfully.')

    data_df.columns = [str(col).strip().lower() for col in data_df.columns]
    data_df = validate.col_header(data_df, cols.WAVELENGTH, validate.WAVELENGTH_ALIASES)
    data_df = validate.col_header(data_df, cols.DARK, validate.DARK_ALIASES)
    data_df, irr_wls = validate.irr_header(data_df)

    validate.violin_irrs(irr_wls, violin_irrs)

    ordered_cols = [cols.WAVELENGTH, cols.DARK] + irr_wls
    data_df = data_df[ordered_cols].copy()

    for col in ordered_cols:
        data_df[col] = pd.to_numeric(data_df[col], errors='coerce')

    validate.numeric_input(data_df, irr_wls)

    wl_data = data_df[cols.WAVELENGTH]

    validate.irradiation_wavelength_range(wl_data, irr_wls)

    if singlepair:
        validate.singlepair(irr_wls, singlepair)
    
    validate.range_overlap_wl(
        wl_data,
        darkmax_range,
        'darkmax_range',
    )

    validate.range_overlap_wl(
        wl_data,
        meta_range,
        'meta_range',
    )

    validate.range_overlap_wl(
        wl_data,
        uvvis_range,
        'uvvis_range',
    )


    data_df = data_df.sort_values(cols.WAVELENGTH, ascending=False).reset_index(drop=True)

    outlog.data_overview(data_df, irr_wls, uvvis_range, darkmax_range)

    return data_df, irr_wls